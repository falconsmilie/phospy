from __future__ import annotations

import socket
from collections.abc import Mapping
from dataclasses import replace

import pandas as pd
import pytest

from phospy.api import (
    EnrichmentConfig,
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
    ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_DROP,
    SUPPORTED_MULTIPLE_TESTING_CORRECTIONS,
)
from phospy.science.enrichment.ora import OraEngine
from phospy.science.statistics.multiple_testing import (
    run as run_multiple_testing_correction,
)
from phospy.workflows import EnrichmentWorkflow
from phospy.workflows.enrichment.executor import EnrichmentWorkflowExecutor
from phospy.workflows.enrichment.interpreter import EnrichmentWorkflowInterpreter
from phospy.workflows.enrichment.validator import EnrichmentWorkflowValidator


def _gene_collection() -> GeneSetCollection:
    return GeneSetCollection(
        sets={
            "KINASE_RESPONSE": ("AKT1", "MAPK1", "MTOR"),
            "CELL_CYCLE": ("CDK1", "CDK2", "MAPK1"),
            "EMPTY_AFTER_BACKGROUND": ("OUTSIDE_A", "OUTSIDE_B"),
        },
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        term_names={
            "KINASE_RESPONSE": "Kinase response",
            "CELL_CYCLE": "Cell cycle",
            "EMPTY_AFTER_BACKGROUND": "Empty after background",
        },
        source_name="unit_test",
        source_version="2026.06",
    )


def _ptm_collection() -> PtmSetCollection:
    return PtmSetCollection(
        sets={
            "MOTIF_A": ("rat|P12345|S10", "rat|P12345|T20"),
            "MOTIF_B": ("rat|P12345|Y30",),
        },
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
        source_name="unit_test",
    )


def test_enrichment_workflow_happy_path_with_gene_symbols() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(),
    )

    result = EnrichmentWorkflow().run(request)

    assert isinstance(result, EnrichmentWorkflowResult)
    assert result.identifier_kind == ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
    assert result.method_metadata["method"] == "over_representation"
    assert result.background_summary["source"] == "explicit"
    assert result.background_summary["universe_size"] == 5
    assert result.set_collection_summary["collection_kind"] == "gene_set"
    foreground_background = result.diagnostics["foreground_background"]
    assert foreground_background["identifier_kind"] == (
        ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
    )
    assert foreground_background["foreground_size_before_intersection"] == 2
    assert foreground_background["background_size"] == 5
    assert (
        foreground_background["usable_foreground_size_after_background_intersection"]
        == 2
    )
    assert foreground_background["foreground_identifiers_missing_from_background"] == ()
    assert foreground_background["tested_set_count"] == 3
    assert foreground_background["dropped_set_count"] == 0
    assert foreground_background["set_identifiers_missing_from_background"] == (
        "OUTSIDE_A",
        "OUTSIDE_B",
    )
    assert result.provenance is not None
    assert tuple(result.table["term_id"])[:1] == ("KINASE_RESPONSE",)


def test_enrichment_selected_outside_background_defaults_to_error() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "UNKNOWN"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="selected_identifiers.*selected_outside_background_policy='error'",
    ):
        EnrichmentWorkflow().run(request)


def test_enrichment_explicit_drop_records_universe_policy_diagnostics() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "UNKNOWN", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(
            selected_outside_background_policy=(
                ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_DROP
            )
        ),
    )

    result = EnrichmentWorkflow().run(request)

    foreground_background = result.diagnostics["foreground_background"]
    assert foreground_background["selected_outside_background_policy"] == "drop"
    assert foreground_background["foreground_identifiers_missing_from_background"] == (
        "UNKNOWN",
    )
    assert foreground_background["retained_foreground_fraction"] == pytest.approx(2 / 3)
    assert result.provenance is not None
    universe_policy = result.provenance.workflow_parameters["universe_policy"]
    assert universe_policy["selected_identifiers_outside_background_count"] == 1
    assert universe_policy["selected_identifiers_outside_background"] == ["UNKNOWN"]
    assert universe_policy["set_identifiers_outside_background"] == [
        "OUTSIDE_A",
        "OUTSIDE_B",
    ]


def test_enrichment_namespace_mismatch_is_not_silently_dropped() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("rat|P12345|S10",),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(),
    )

    with pytest.raises(WorkflowValidationError, match="namespace mismatch"):
        EnrichmentWorkflow().run(request)


def test_enrichment_set_member_outside_background_policy_is_configurable() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(set_member_outside_background_policy="error"),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="set_collection identifiers.*set_member_outside_background_policy",
    ):
        EnrichmentWorkflow().run(request)


def test_enrichment_retained_foreground_fraction_threshold_is_enforced() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1", "UNKNOWN"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(
            selected_outside_background_policy=(
                ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_DROP
            ),
            minimum_retained_foreground_fraction=0.75,
        ),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="retained foreground fraction.*minimum_retained_foreground_fraction",
    ):
        EnrichmentWorkflow().run(request)


def test_enrichment_set_outside_background_diagnostics_are_order_deterministic() -> (
    None
):
    collection_a = GeneSetCollection(
        sets={
            "B_TERM": ("OUTSIDE_B", "AKT1"),
            "A_TERM": ("MAPK1", "OUTSIDE_A"),
        },
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        source_name="unit_test",
    )
    collection_b = GeneSetCollection(
        sets={
            "A_TERM": ("OUTSIDE_A", "MAPK1"),
            "B_TERM": ("AKT1", "OUTSIDE_B"),
        },
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        source_name="unit_test",
    )

    def _run(collection: GeneSetCollection) -> EnrichmentWorkflowResult:
        return EnrichmentWorkflow().run(
            EnrichmentWorkflowRequest(
                identifier_column="gene_symbol",
                identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
                set_collection=collection,
                selected_identifiers=("AKT1", "MAPK1"),
                background_universe=("AKT1", "MAPK1", "MTOR"),
                config=EnrichmentConfig(),
            )
        )

    result_a = _run(collection_a)
    result_b = _run(collection_b)

    assert result_a.diagnostics["foreground_background"][
        "set_identifiers_missing_from_background"
    ] == ("OUTSIDE_A", "OUTSIDE_B")
    assert result_b.diagnostics["foreground_background"][
        "set_identifiers_missing_from_background"
    ] == ("OUTSIDE_A", "OUTSIDE_B")
    assert result_a.provenance is not None
    assert result_b.provenance is not None
    assert (
        result_a.provenance.workflow_parameters["universe_policy"]
        == (result_b.provenance.workflow_parameters["universe_policy"])
    )


def test_enrichment_workflow_provenance_records_method_and_limitations() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(
            multiple_testing_correction="benjamini_hochberg",
        ),
    )

    result = EnrichmentWorkflow().run(request)

    provenance = result.provenance
    assert provenance is not None
    parameters = provenance.workflow_parameters
    assert parameters["method"] == "over_representation"
    assert parameters["identifier_kind"] == ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
    assert parameters["background_universe_size"] == 5
    assert parameters["selected_identifier_count"] == 2
    assert parameters["multiple_testing_correction"] == "benjamini_hochberg"
    assert parameters["multiple_testing_method"] == "benjamini_hochberg"
    assert parameters["number_of_tests"] == len(result.records)
    assert parameters["correction_owner"] == "ora_engine"
    assert parameters["offline_no_online_resource_policy"] == (
        "offline_user_supplied_collections_only"
    )
    assert parameters["online_resources_used"] is False

    set_collection = parameters["set_collection"]
    assert isinstance(set_collection, dict)
    assert set_collection["source_name"] == "unit_test"
    assert set_collection["source_version"] == "2026.06"
    limitations = tuple(parameters["limitations"])
    assert "offline over-representation analysis only" in limitations
    assert any("GO, KEGG, Reactome, and PTM-SEA" in item for item in limitations)


def test_enrichment_workflow_provenance_records_identifier_and_set_attrition() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={
                "TESTED": ("AKT1", "MAPK1"),
                "TOO_SMALL": ("OUTSIDE_A",),
                "TOO_LARGE": ("AKT1", "MAPK1", "MTOR", "CDK1"),
            },
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            source_name="unit_test",
        ),
        selected_identifiers=("AKT1", "AKT1", "OUTSIDE_SELECTED", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK1"),
        config=EnrichmentConfig(
            min_set_size=1,
            max_set_size=3,
            selected_outside_background_policy=(
                ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_DROP
            ),
        ),
    )

    result = EnrichmentWorkflow().run(request)

    assert result.provenance is not None
    parameters = result.provenance.workflow_parameters
    metrics = parameters["row_attrition_metrics"]
    assert metrics["selected_identifiers_provided"] == 4
    assert metrics["selected_identifiers_prepared"] == 3
    assert metrics["selected_identifiers_retained_in_universe"] == 2
    assert metrics["selected_identifiers_dropped_outside_universe"] == 1
    assert metrics["background_identifiers_provided"] == 5
    assert metrics["background_identifiers_retained_in_universe"] == 4
    assert metrics["sets_provided"] == 3
    assert metrics["sets_tested"] == 1
    assert metrics["sets_skipped_due_to_min_max_size"] == 2
    assert metrics["sets_skipped_due_to_min_size"] == 1
    assert metrics["sets_skipped_due_to_max_size"] == 1

    row_attrition = parameters["row_attrition"]
    selected = row_attrition["selected_identifiers"]
    assert selected["input_rows"] == 4
    assert selected["final_rows"] == 2
    assert [record["reason"] for record in selected["records"]] == [
        "selected_identifiers_dropped_before_universe_intersection",
        "selected_identifiers_not_retained_in_universe",
    ]
    background = row_attrition["background_identifiers"]
    assert background["input_rows"] == 5
    assert background["final_rows"] == 4
    sets = row_attrition["sets"]
    assert sets["input_rows"] == 3
    assert sets["final_rows"] == 1
    assert sets["records"][0]["examples"] == ["TOO_SMALL", "TOO_LARGE"]


def test_enrichment_workflow_happy_path_with_site_keys() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="site_key",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
        set_collection=_ptm_collection(),
        input_table=pd.DataFrame({"site_key": ["rat|P12345|S10", " rat|P12345|T20 "]}),
        background_universe=(
            "rat|P12345|S10",
            "rat|P12345|T20",
            "rat|P12345|Y30",
        ),
    )

    result = EnrichmentWorkflow().run(request)

    assert result.identifier_kind == ENRICHMENT_IDENTIFIER_KIND_SITE_KEY
    assert result.set_collection_summary["collection_kind"] == "ptm_set"
    assert result.background_summary["selected_identifier_source"] == "input_table"
    assert result.diagnostics["foreground_background"]["identifier_kind"] == (
        ENRICHMENT_IDENTIFIER_KIND_SITE_KEY
    )
    assert result.records[0].overlap_identifiers == (
        "rat|P12345|S10",
        "rat|P12345|T20",
    )


def test_enrichment_empty_foreground_after_explicit_drop_is_rejected() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("UNKNOWN_A", "UNKNOWN_B"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(
            selected_outside_background_policy=(
                ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_DROP
            )
        ),
    )

    with pytest.raises(
        WorkflowValidationError,
        match="empty after background filtering",
    ):
        EnrichmentWorkflow().run(request)


def test_enrichment_mixed_matched_and_unmatched_foreground_keeps_ora_statistics() -> (
    None
):
    matched_request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(),
    )
    mixed_request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "UNKNOWN", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(
            selected_outside_background_policy=(
                ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_DROP
            )
        ),
    )

    matched_result = EnrichmentWorkflow().run(matched_request)
    mixed_result = EnrichmentWorkflow().run(mixed_request)

    pd.testing.assert_frame_equal(mixed_result.table, matched_result.table)
    foreground_background = mixed_result.diagnostics["foreground_background"]
    assert foreground_background["foreground_size_before_intersection"] == 3
    assert (
        foreground_background["usable_foreground_size_after_background_intersection"]
        == 2
    )
    assert foreground_background["foreground_identifiers_missing_from_background"] == (
        "UNKNOWN",
    )


def test_enrichment_workflow_stage_collaboration() -> None:
    calls: list[tuple[str, object]] = []
    request = object()
    validated = object()
    interpreted = object()
    expected = object()

    class Validator:
        def run(self, value: object) -> object:
            calls.append(("validator", value))
            return validated

    class Interpreter:
        def run(self, value: object) -> object:
            calls.append(("interpreter", value))
            return interpreted

    class Executor:
        def run(self, value: object) -> object:
            calls.append(("executor", value))
            return expected

    result = EnrichmentWorkflow._with_components(
        validator=Validator(),
        interpreter=Interpreter(),
        executor=Executor(),
    ).run(request)  # type: ignore[arg-type]

    assert result is expected
    assert calls == [
        ("validator", request),
        ("interpreter", validated),
        ("executor", interpreted),
    ]


def test_enrichment_result_table_shape() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
    )

    table = EnrichmentWorkflow().run(request).table

    assert tuple(table.columns) == (
        "term_id",
        "term_name",
        "collection_kind",
        "identifier_kind",
        "input_overlap_count",
        "background_overlap_count",
        "set_size",
        "overlap_identifiers",
        "p_value",
        "adjusted_p_value",
        "correction_method",
        "enrichment_ratio",
    )
    assert table.shape[0] == 3


def test_enrichment_workflow_applies_multiple_testing_correction() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
    )

    result = EnrichmentWorkflow().run(request)

    assert all(record.adjusted_p_value is not None for record in result.records)
    assert all(record.p_value is not None for record in result.records)
    assert all(
        record.correction_method == "benjamini_hochberg" for record in result.records
    )
    assert result.diagnostics["multiple_testing_correction"] == {
        "method": "benjamini_hochberg",
        "applied": True,
        "tested_record_count": 3,
    }


def test_enrichment_workflow_reports_none_correction_not_applied() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(multiple_testing_correction="none"),
    )

    result = EnrichmentWorkflow().run(request)

    assert result.diagnostics["multiple_testing_correction"] == {
        "method": "none",
        "applied": False,
        "tested_record_count": 3,
    }
    assert tuple(record.adjusted_p_value for record in result.records) == (
        pytest.approx(tuple(record.p_value for record in result.records))
    )
    assert all(record.correction_method == "none" for record in result.records)


@pytest.mark.parametrize("method", SUPPORTED_MULTIPLE_TESTING_CORRECTIONS)
def test_enrichment_workflow_supports_configured_correction_methods(
    method: MultipleTestingCorrection,
) -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(multiple_testing_correction=method),
    )

    result = EnrichmentWorkflow().run(request)

    raw_p_values = tuple(record.p_value for record in result.records)
    expected_adjusted = run_multiple_testing_correction(
        raw_p_values,
        method=method,
    )
    assert tuple(record.adjusted_p_value for record in result.records) == (
        pytest.approx(expected_adjusted)
    )
    assert all(record.correction_method == method for record in result.records)
    assert result.method_metadata["multiple_testing_correction"] == method
    correction_diagnostics = result.diagnostics["multiple_testing_correction"]
    expected_applied = method != "none"
    assert isinstance(correction_diagnostics, Mapping)
    assert correction_diagnostics["method"] == method
    assert correction_diagnostics["applied"] is expected_applied
    assert correction_diagnostics["tested_record_count"] == len(result.records)
    assert result.provenance is not None
    parameters = result.provenance.workflow_parameters
    assert parameters["multiple_testing_correction"] == method
    assert parameters["multiple_testing_method"] == method
    assert parameters["number_of_tests"] == len(result.records)
    assert parameters["correction_owner"] == "ora_engine"


def test_enrichment_workflow_empty_and_no_hit_sets_are_reported() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1",),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
    )

    result = EnrichmentWorkflow().run(request)
    rows = {record.term_id: record for record in result.records}

    assert rows["EMPTY_AFTER_BACKGROUND"].background_overlap_count == 0
    assert rows["EMPTY_AFTER_BACKGROUND"].input_overlap_count == 0
    assert rows["EMPTY_AFTER_BACKGROUND"].p_value == pytest.approx(1.0)
    assert result.set_collection_summary["empty_after_background_count"] == 1
    assert "CELL_CYCLE" in result.diagnostics["ora"]["no_hit_set_ids"]


def test_enrichment_executor_uses_ora_adjusted_p_values_without_overwrite() -> None:
    validated = EnrichmentWorkflowValidator().run(
        EnrichmentWorkflowRequest(
            identifier_column="gene_symbol",
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            set_collection=_gene_collection(),
            selected_identifiers=("AKT1", "MAPK1"),
            background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        )
    )
    interpreted = EnrichmentWorkflowInterpreter().run(validated)
    calls: list[str] = []
    engine_adjusted_values = (0.11, 0.22, 0.33)

    class Engine:
        def run(self, **kwargs: object) -> object:
            calls.append("ora")
            assert kwargs["selected_identifiers"] == ("AKT1", "MAPK1")
            assert kwargs["config"].multiple_testing_correction == "benjamini_hochberg"
            ora_result = OraEngine().run(**kwargs)  # type: ignore[arg-type]
            adjusted_records = tuple(
                replace(
                    record,
                    adjusted_p_value=adjusted_p_value,
                )
                for record, adjusted_p_value in zip(
                    ora_result.records,
                    engine_adjusted_values,
                    strict=True,
                )
            )
            return replace(ora_result, records=adjusted_records)

    result = EnrichmentWorkflowExecutor(
        ora_engine=Engine(),  # type: ignore[arg-type]
    ).run(interpreted)

    assert calls == ["ora"]
    assert (
        tuple(record.adjusted_p_value for record in result.records)
        == engine_adjusted_values
    )
    assert all(record.p_value is not None for record in result.records)


def test_enrichment_workflow_has_no_internet_dependency(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_network(*args: object, **kwargs: object) -> object:
        raise AssertionError("enrichment workflow must run offline")

    monkeypatch.setattr(socket, "create_connection", fail_network)
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
    )

    result = EnrichmentWorkflow().run(request)

    assert result.table.shape[0] == 3


def test_enrichment_set_size_filter_excludes_sets_below_minimum() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={
                "SMALL": ("AKT1",),
                "PASS": ("AKT1", "MAPK1"),
            },
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            source_name="unit_test",
        ),
        selected_identifiers=("AKT1",),
        background_universe=("AKT1", "MAPK1", "MTOR"),
        config=EnrichmentConfig(min_set_size=2),
    )

    result = EnrichmentWorkflow().run(request)

    assert tuple(record.term_id for record in result.records) == ("PASS",)
    assert result.set_collection_summary["dropped_set_count"] == 1
    assert result.set_collection_summary["dropped_set_ids"] == ("SMALL",)
    diagnostics = result.diagnostics["set_size_filter"]
    assert isinstance(diagnostics, Mapping)
    assert diagnostics["dropped_set_reason_counts"] == {
        "below_min_set_size": 1,
        "above_max_set_size": 0,
    }
    assert diagnostics["dropped_sets"] == (
        {
            "set_id": "SMALL",
            "term_name": "SMALL",
            "reason": "below_min_set_size",
            "raw_set_size": 1,
            "background_overlap_count": 1,
            "identifiers_outside_background_count": 0,
        },
    )
    foreground_background = result.diagnostics["foreground_background"]
    assert foreground_background["tested_set_count"] == 1
    assert foreground_background["dropped_set_count"] == 1


def test_enrichment_set_size_filter_excludes_sets_above_maximum() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={
                "LARGE": ("AKT1", "MAPK1", "MTOR"),
                "PASS": ("AKT1", "MAPK1"),
            },
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            source_name="unit_test",
        ),
        selected_identifiers=("AKT1",),
        background_universe=("AKT1", "MAPK1", "MTOR"),
        config=EnrichmentConfig(max_set_size=2),
    )

    result = EnrichmentWorkflow().run(request)

    assert tuple(record.term_id for record in result.records) == ("PASS",)
    assert result.set_collection_summary["dropped_set_count"] == 1
    assert result.set_collection_summary["dropped_set_reason_counts"] == {
        "below_min_set_size": 0,
        "above_max_set_size": 1,
    }


def test_enrichment_set_size_filter_uses_background_intersection_for_passing_set() -> (
    None
):
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={"PASS_AFTER_BACKGROUND": ("AKT1", "MAPK1", "OUTSIDE")},
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            source_name="unit_test",
        ),
        selected_identifiers=("AKT1",),
        background_universe=("AKT1", "MAPK1", "MTOR"),
        config=EnrichmentConfig(min_set_size=2),
    )

    result = EnrichmentWorkflow().run(request)

    assert tuple(record.term_id for record in result.records) == (
        "PASS_AFTER_BACKGROUND",
    )
    assert result.records[0].set_size == 3
    assert result.records[0].background_overlap_count == 2
    assert result.set_collection_summary["identifiers_outside_background_count"] == 1
    assert result.diagnostics["set_size_filter"]["dropped_set_count"] == 0


def test_enrichment_set_size_filter_can_fail_only_after_background_intersection() -> (
    None
):
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={
                "FAIL_AFTER_BACKGROUND": ("AKT1", "OUTSIDE_A", "OUTSIDE_B"),
                "PASS": ("AKT1", "MAPK1"),
            },
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            source_name="unit_test",
        ),
        selected_identifiers=("AKT1",),
        background_universe=("AKT1", "MAPK1", "MTOR"),
        config=EnrichmentConfig(min_set_size=2),
    )

    result = EnrichmentWorkflow().run(request)

    assert tuple(record.term_id for record in result.records) == ("PASS",)
    dropped_sets = result.diagnostics["set_size_filter"]["dropped_sets"]
    assert dropped_sets == (
        {
            "set_id": "FAIL_AFTER_BACKGROUND",
            "term_name": "FAIL_AFTER_BACKGROUND",
            "reason": "below_min_set_size",
            "raw_set_size": 3,
            "background_overlap_count": 1,
            "identifiers_outside_background_count": 2,
        },
    )


def test_enrichment_set_size_filter_reports_all_sets_dropped() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={
                "SMALL": ("AKT1",),
                "OUTSIDE_ONLY": ("OUTSIDE",),
            },
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            source_name="unit_test",
        ),
        selected_identifiers=("AKT1",),
        background_universe=("AKT1", "MAPK1", "MTOR"),
        config=EnrichmentConfig(min_set_size=2),
    )

    result = EnrichmentWorkflow().run(request)

    assert result.records == ()
    assert result.table.empty
    assert result.unmatched_identifiers == ("AKT1",)
    assert result.set_collection_summary["tested_set_count"] == 0
    assert result.set_collection_summary["dropped_set_count"] == 2
    assert result.diagnostics["ora"]["record_count"] == 0
    assert result.diagnostics["multiple_testing_correction"] == {
        "method": "benjamini_hochberg",
        "applied": False,
        "tested_record_count": 0,
    }


def test_enrichment_set_size_filter_adjusts_p_values_over_tested_sets_only() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={
                "TESTED": ("A", "B", "D", "E"),
                "DROPPED_SMALL": ("A",),
            },
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            source_name="unit_test",
        ),
        selected_identifiers=("A", "B", "C"),
        background_universe=("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
        config=EnrichmentConfig(min_set_size=2),
    )

    result = EnrichmentWorkflow().run(request)

    assert tuple(record.term_id for record in result.records) == ("TESTED",)
    assert result.records[0].p_value == pytest.approx(1.0 / 3.0)
    assert result.records[0].adjusted_p_value == pytest.approx(1.0 / 3.0)
    assert result.diagnostics["multiple_testing_correction"]["tested_record_count"] == 1


def test_enrichment_correction_denominator_ignores_dropped_sets() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={
                "TESTED_STRONG": ("A", "B", "C"),
                "TESTED_WEAK": ("A", "D", "E", "F"),
                "DROPPED_SMALL": ("A",),
            },
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            source_name="unit_test",
        ),
        selected_identifiers=("A", "B", "C"),
        background_universe=("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"),
        config=EnrichmentConfig(
            min_set_size=2,
            multiple_testing_correction="bonferroni",
        ),
    )

    result = EnrichmentWorkflow().run(request)

    rows = {record.term_id: record for record in result.records}
    assert set(rows) == {"TESTED_STRONG", "TESTED_WEAK"}
    assert rows["TESTED_STRONG"].p_value == pytest.approx(1.0 / 120.0)
    assert rows["TESTED_STRONG"].adjusted_p_value == pytest.approx(1.0 / 60.0)
    assert result.diagnostics["multiple_testing_correction"]["tested_record_count"] == 2
    assert result.set_collection_summary["dropped_set_ids"] == ("DROPPED_SMALL",)


def test_enrichment_set_size_filter_defaults_leave_results_unchanged() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(),
    )
    explicit_default_request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR", "CDK1", "CDK2"),
        config=EnrichmentConfig(min_set_size=None, max_set_size=None),
    )

    result = EnrichmentWorkflow().run(request)
    explicit_default_result = EnrichmentWorkflow().run(explicit_default_request)

    pd.testing.assert_frame_equal(result.table, explicit_default_result.table)
    assert "set_size_filter" not in result.diagnostics
    assert "set_size_filter" not in explicit_default_result.diagnostics

from __future__ import annotations

import socket

import pandas as pd
import pytest

from phospy.api import (
    ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
    EnrichmentConfig,
    EnrichmentWorkflowRequest,
    EnrichmentWorkflowResult,
    GeneSetCollection,
    PtmSetCollection,
)
from phospy.science.enrichment.ora import OraEngine
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
    assert result.provenance is not None
    assert tuple(result.table["term_id"])[:1] == ("KINASE_RESPONSE",)


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
    assert result.records[0].overlap_identifiers == (
        "rat|P12345|S10",
        "rat|P12345|T20",
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

    result = EnrichmentWorkflow(
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
    assert all(
        record.correction_method == "benjamini_hochberg" for record in result.records
    )
    assert result.diagnostics["multiple_testing_correction"] == {
        "method": "benjamini_hochberg",
        "applied": True,
        "tested_record_count": 3,
    }


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


def test_enrichment_executor_calls_ora_and_correction_helper() -> None:
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

    class Engine:
        def run(self, **kwargs: object) -> object:
            calls.append("ora")
            assert kwargs["selected_identifiers"] == ("AKT1", "MAPK1")
            return OraEngine().run(**kwargs)  # type: ignore[arg-type]

    def correction_runner(
        p_values: tuple[float | None, ...],
        *,
        method: str,
    ) -> tuple[float | None, ...]:
        calls.append("correction")
        assert method == "benjamini_hochberg"
        return tuple(0.5 for _ in p_values)

    result = EnrichmentWorkflowExecutor(
        ora_engine=Engine(),  # type: ignore[arg-type]
        correction_runner=correction_runner,  # type: ignore[arg-type]
    ).run(interpreted)

    assert calls == ["ora", "correction"]
    assert tuple(record.adjusted_p_value for record in result.records) == (
        0.5,
        0.5,
        0.5,
    )


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

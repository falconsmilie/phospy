from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    MULTIPLE_TESTING_CORRECTION_NONE,
    EnrichmentConfig,
    EnrichmentWorkflowRequest,
    GeneSetCollection,
)
from phospy.errors import WorkflowBoundaryError
from phospy.science.enrichment.ora import OraResult, OraResultRecord
from phospy.workflows.enrichment.executor import EnrichmentWorkflowExecutor
from phospy.workflows.enrichment.interpreter import EnrichmentWorkflowInterpreter
from phospy.workflows.enrichment.models import (
    InterpretedEnrichmentWorkflowRequest,
    ValidatedEnrichmentWorkflowRequest,
)
from phospy.workflows.enrichment.validator import EnrichmentWorkflowValidator
from phospy.workflows.enrichment.workflow import EnrichmentWorkflow


def _gene_collection() -> GeneSetCollection:
    return GeneSetCollection(
        sets={
            "MAPK_PATHWAY": ("AKT1", "MAPK1"),
            "MTOR_SIGNALING": ("MTOR",),
        },
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        term_names={
            "MAPK_PATHWAY": "MAPK pathway",
            "MTOR_SIGNALING": "MTOR signaling",
        },
        source_name="unit_test",
    )


def _valid_gene_request() -> EnrichmentWorkflowRequest:
    return EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", "MAPK1"),
        background_universe=("AKT1", "MAPK1", "MTOR"),
        config=EnrichmentConfig(),
    )


def test_enrichment_workflow_calls_validator_interpreter_executor_in_order() -> None:
    events: list[str] = []
    validated = object()
    interpreted = object()
    expected_result = object()

    class _Validator:
        def run(self, request: object) -> object:
            events.append("validator")
            return validated

    class _Interpreter:
        def run(self, request: object) -> object:
            events.append("interpreter")
            assert request is validated
            return interpreted

    class _Executor:
        def run(self, request: object) -> object:
            events.append("executor")
            assert request is interpreted
            return expected_result

    result = EnrichmentWorkflow(
        validator=_Validator(),  # type: ignore[arg-type]
        interpreter=_Interpreter(),  # type: ignore[arg-type]
        executor=_Executor(),  # type: ignore[arg-type]
    ).run(object())  # type: ignore[arg-type]

    assert events == ["validator", "interpreter", "executor"]
    assert result is expected_result


def test_enrichment_stage_components_expose_run() -> None:
    for stage_type in (
        EnrichmentWorkflow,
        EnrichmentWorkflowValidator,
        EnrichmentWorkflowInterpreter,
        EnrichmentWorkflowExecutor,
    ):
        assert callable(getattr(stage_type, "run", None))


def test_enrichment_workflow_validator_does_not_run_ora(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.science.enrichment.ora as ora

    def fail_if_called(*args: object, **kwargs: object) -> object:
        raise AssertionError("ORA must not run during enrichment validation")

    monkeypatch.setattr(ora, "run", fail_if_called)
    monkeypatch.setattr(ora.OraEngine, "run", fail_if_called)

    validated = EnrichmentWorkflowValidator().run(_valid_gene_request())

    assert isinstance(validated, ValidatedEnrichmentWorkflowRequest)
    assert validated.selected_identifiers == ("AKT1", "MAPK1")
    assert not hasattr(validated, "method_config")
    assert not hasattr(validated, "diagnostics")


def test_enrichment_interpreter_resolves_set_and_background_execution_inputs() -> None:
    validated = EnrichmentWorkflowValidator().run(_valid_gene_request())

    interpreted = EnrichmentWorkflowInterpreter().run(validated)

    assert isinstance(interpreted, InterpretedEnrichmentWorkflowRequest)
    assert interpreted.selected_identifiers == validated.selected_identifiers
    assert interpreted.background_universe == validated.background_universe
    assert interpreted.set_collection is validated.set_collection
    assert interpreted.method_config.statistical_test == "hypergeometric"
    assert interpreted.method_config.selected_outside_background_policy == "drop"
    assert interpreted.method_config.set_outside_background_policy == "drop"
    assert interpreted.identifier_semantics.analysis_level == "gene"
    assert interpreted.background_summary == {
        "source": "explicit",
        "universe_size": 3,
        "selected_identifier_count": 2,
        "selected_identifier_source": "selected_identifiers",
    }
    assert interpreted.set_collection_summary["set_count"] == 2
    assert interpreted.set_collection_summary["distinct_member_count"] == 3
    assert interpreted.diagnostics["interpreter"] == {
        "selected_identifiers_prepared": 2,
        "background_universe_prepared": 3,
        "set_count_prepared": 2,
    }


def test_enrichment_executor_accepts_only_interpreted_requests() -> None:
    with pytest.raises(
        WorkflowBoundaryError,
        match="enrichment.executor.interpreted_request_type",
    ):
        EnrichmentWorkflowExecutor().run(_valid_gene_request())  # type: ignore[arg-type]


def test_enrichment_executor_consumes_interpreted_inputs_without_inference() -> None:
    interpreted = EnrichmentWorkflowInterpreter().run(
        EnrichmentWorkflowValidator().run(_valid_gene_request())
    )
    captured: dict[str, object] = {}

    class _OraEngineSpy:
        def run(
            self,
            *,
            selected_identifiers,
            background_universe,
            enrichment_sets,
            config,
        ) -> OraResult:
            captured["selected_identifiers"] = tuple(selected_identifiers)
            captured["background_universe"] = tuple(background_universe)
            captured["enrichment_sets"] = enrichment_sets
            captured["config"] = config
            return OraResult(
                method="over_representation",
                config=config,
                background_size=len(tuple(background_universe)),
                selected_size=len(tuple(selected_identifiers)),
                selected_identifiers=tuple(selected_identifiers),
                dropped_selected_identifiers=(),
                records=(
                    OraResultRecord(
                        set_id="MAPK_PATHWAY",
                        name="MAPK pathway",
                        collection_kind=interpreted.set_collection.collection_kind,
                        identifier_kind=interpreted.identifier_semantics.identifier_kind,
                        background_size=len(interpreted.background_universe),
                        selected_size=len(interpreted.selected_identifiers),
                        raw_set_size=2,
                        set_size=2,
                        overlap_size=2,
                        overlap_identifiers=("AKT1", "MAPK1"),
                        p_value=0.25,
                        enrichment_ratio=1.5,
                        set_identifiers_outside_background_count=0,
                    ),
                ),
            )

    def correction_runner(p_values, *, method):
        captured["correction_p_values"] = tuple(p_values)
        captured["correction_method"] = method
        return (0.5,)

    result = EnrichmentWorkflowExecutor(
        ora_engine=_OraEngineSpy(),
        correction_runner=correction_runner,
    ).run(interpreted)

    assert captured["selected_identifiers"] == interpreted.selected_identifiers
    assert captured["background_universe"] == interpreted.background_universe
    assert captured["enrichment_sets"] is interpreted.set_collection
    assert captured["config"].multiple_testing_correction == (
        MULTIPLE_TESTING_CORRECTION_NONE
    )
    assert captured["correction_p_values"] == (0.25,)
    assert captured["correction_method"] == (
        interpreted.method_config.multiple_testing_correction
    )
    assert result.records[0].adjusted_p_value == 0.5
    assert result.background_summary["source"] == "explicit"
    assert isinstance(result.table, pd.DataFrame)

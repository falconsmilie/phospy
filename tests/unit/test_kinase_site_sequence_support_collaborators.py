from __future__ import annotations

import inspect
from typing import cast

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset, KinaseWorkflow
from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.workflows.kinase.contracts import (
    KinaseWorkflowExecutorContract,
    ResolvedKinaseExecutionConfig,
    ResolvedKinaseWorkflowRequest,
)
from phospy.workflows.kinase.interpreter import KinaseWorkflowInterpreter
from phospy.workflows.kinase.scoring_runner import KinaseScoringRunner
from phospy.workflows.kinase.site_sequence_support import (
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    KinaseSiteSequenceSupportBuilder,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def _dataset(
    *,
    sequences_by_site: dict[str, str],
) -> AnalysisReadyPhosphoDataset:
    site_ids = list(sequences_by_site)
    site_index = pd.Index(site_ids, name="site_id")
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {
                "sample_a": [1.0 + float(i) for i, _ in enumerate(site_ids)],
                "sample_b": [2.0 + float(i) for i, _ in enumerate(site_ids)],
            },
            index=site_index,
        ),
        site_metadata=pd.DataFrame(
            {
                "gene_symbol": [site.split(";", 1)[0] for site in site_ids],
                "site": [site.split(";")[1] for site in site_ids],
                "site_sequence": [sequences_by_site[site] for site in site_ids],
            },
            index=site_index,
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references(
    *,
    reference_sequences_by_site: dict[str, str],
) -> ReferenceBundle:
    site_ids = list(reference_sequences_by_site)
    site_index = pd.Index(site_ids, name="site_id")
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6" for _ in site_ids],
                "substrate_site": site_ids,
            }
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    reference_sequences_by_site[site_id] for site_id in site_ids
                ]
            },
            index=site_index,
        ),
    )


class _ExecutorSentinel:
    def __init__(self) -> None:
        self.called = False

    def run(self, request: ResolvedKinaseWorkflowRequest) -> object:
        _ = request
        self.called = True
        raise AssertionError("executor should not run for interpreter conflict errors")


def test_site_sequence_builder_appends_dataset_only_sequences() -> None:
    dataset = _dataset(
        sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA",
            "EXTRA;S1;": "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )

    result = KinaseSiteSequenceSupportBuilder().run(
        dataset=dataset.phospho,
        site_metadata=dataset.site_metadata,
        reference_site_sequences=references.site_sequences,
        conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    )

    assert "EXTRA;S1;" in result.site_sequences.index
    assert (
        result.site_sequences.at["EXTRA;S1;", "site_sequence"]
        == "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA"
    )
    assert result.dataset_sequences_added == 1
    assert result.dataset_reference_conflict_count == 0
    assert result.conflicts == ()


def test_site_sequence_builder_ignores_matching_dataset_reference_sequences() -> None:
    dataset = _dataset(
        sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )

    result = KinaseSiteSequenceSupportBuilder().run(
        dataset=dataset.phospho,
        site_metadata=dataset.site_metadata,
        reference_site_sequences=references.site_sequences,
        conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    )

    assert result.dataset_sequences_added == 0
    assert result.dataset_reference_conflict_count == 0
    assert result.conflicts == ()


def test_interpreter_conflict_error_policy_fails_before_executor_runs() -> None:
    dataset = _dataset(
        sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAATTTTTTTTTTTTTTTTTTTTTTTT",
            "GSK3B;S9;": "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )
    executor = _ExecutorSentinel()

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow(
            interpreter=KinaseWorkflowInterpreter(
                site_sequence_conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR
            ),
            executor=cast(KinaseWorkflowExecutorContract, executor),
        ).run(request)

    error = exc_info.value
    assert error.seam == "kinase.interpreter.site_sequence_conflict"
    diagnostics = cast(list[dict[str, object]], error.details["conflict_diagnostics"])
    assert len(diagnostics) == 1
    assert diagnostics[0]["site_id"] == "MAPK14;Y182;"
    assert diagnostics[0]["dataset_sequence"] == "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA"
    assert diagnostics[0]["reference_sequence"] == "AAAAAAATTTTTTTTTTTTTTTTTTTTTTTT"
    assert executor.called is False


def test_interpreter_prefer_reference_records_conflict_diagnostics_in_provenance() -> (
    None
):
    dataset = _dataset(
        sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAATTTTTTTTTTTTTTTTTTTTTTTT",
            "GSK3B;S9;": "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )
    result = KinaseWorkflow(
        interpreter=KinaseWorkflowInterpreter(
            site_sequence_conflict_policy=(
                KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE
            )
        )
    ).run(
        KinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(
                top_k=2,
                deterministic_max_selected_kinases=2,
                adaptive_ensemble_runs=2,
            ),
            activity_config=None,
        )
    )

    assert result.provenance is not None
    scoring_diagnostics = cast(
        dict[str, object],
        result.provenance.workflow_parameters["scoring_diagnostics"],
    )
    site_sequence_merge = cast(
        dict[str, object], scoring_diagnostics["site_sequence_merge"]
    )
    assert site_sequence_merge["conflict_policy"] == "prefer_reference"
    conflict_rows = cast(
        list[dict[str, object]], site_sequence_merge["conflict_diagnostics"]
    )
    assert len(conflict_rows) == 1
    assert conflict_rows[0]["selected_sequence_source"] == "reference"
    assert conflict_rows[0]["selected_sequence"] == "AAAAAAATTTTTTTTTTTTTTTTTTTTTTTT"


def test_interpreter_prefer_dataset_selects_dataset_sequence_and_contract_accepts_it() -> (
    None
):
    dataset = _dataset(
        sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAATTTTTTTTTTTTTTTTTTTTTTTT",
            "GSK3B;S9;": "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )
    request = KinaseWorkflowRequest(
        dataset=dataset,
        references=references,
        scoring_config=KinaseScoringConfig(min_substrates=2),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
        ),
        activity_config=None,
    )

    interpreted = KinaseWorkflowInterpreter(
        site_sequence_conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET
    ).run(request)

    assert interpreted.site_sequences.at["MAPK14;Y182;", "site_sequence"] == (
        "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA"
    )
    diagnostics = interpreted.site_sequence_merge_diagnostics
    assert diagnostics["conflict_policy"] == "prefer_dataset"
    conflict_rows = cast(list[dict[str, object]], diagnostics["conflict_diagnostics"])
    assert conflict_rows[0]["selected_sequence_source"] == "dataset"


def test_boundary_scoring_runner_contains_no_site_sequence_merge_conflict_logic() -> (
    None
):
    scoring_runner_source = inspect.getsource(KinaseScoringRunner)
    assert "conflict_policy" not in scoring_runner_source
    assert "conflict_diagnostics" not in scoring_runner_source
    assert "dataset_reference_conflict_count" not in scoring_runner_source


def test_boundary_contracts_reject_unnormalised_site_identifiers() -> None:
    dataset = _dataset(
        sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAYAAAAAAAAAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAASAAAAAAAAAAAAAAAAAAAAAAA",
        }
    )
    unnormalised_sequences = references.site_sequences.copy(deep=True)
    unnormalised_sequences.index = pd.Index(
        [" mapk14 ; y182 ", "GSK3B;S9;"],
        name=references.site_sequences.index.name,
    )

    with pytest.raises(
        WorkflowBoundaryError, match="kinase.contracts.site_sequence_schema"
    ):
        ResolvedKinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            kinase_substrate_map=references.kinase_substrate_map,
            site_sequences=unnormalised_sequences,
            scoring_site_index=dataset.phospho.index.copy(),
            activity_phospho_matrix=dataset.phospho.copy(deep=True),
            execution_config=ResolvedKinaseExecutionConfig(
                scoring_min_substrates=2,
                include_diagnostic_scoring_tables=False,
                profile_missing_value_strategy="strict",
                prediction_top_k=2,
                prediction_deterministic_max_selected_kinases=2,
                prediction_adaptive_ensemble_runs=2,
                prediction_mode="deterministic_ranking",
                prediction_adaptive_policy="stable",
                prediction_n_iterations=5,
                prediction_random_state=None,
            ),
        )

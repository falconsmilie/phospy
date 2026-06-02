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
    KinaseSiteSequenceSupportBuilder,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_from_display_id,
    site_key_index_from_display_ids,
)


def _dataset(
    *,
    sequences_by_site: dict[str, str],
) -> AnalysisReadyPhosphoDataset:
    site_ids = list(sequences_by_site)
    site_index = site_key_index_from_display_ids(site_ids)
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
                "site_key": site_index.astype(str).tolist(),
                "display_id": site_ids,
                **site_key_context_columns(site_index),
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


def _site_identity_map(dataset: AnalysisReadyPhosphoDataset) -> pd.DataFrame:
    metadata = dataset.site_metadata.reindex(dataset.phospho.index)
    site_keys = dataset.phospho.index.astype(str).tolist()
    return pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": metadata.loc[:, "display_id"].astype(str).tolist(),
        },
        index=pd.Index(site_keys, name=dataset.phospho.index.name),
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
            "MAPK14;Y182;": "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "EXTRA;S1;": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
        }
    )

    result = KinaseSiteSequenceSupportBuilder().run(
        dataset=dataset.phospho,
        site_metadata=dataset.site_metadata,
        reference_site_sequences=references.site_sequences,
        conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    )

    extra_site_key = site_key_from_display_id("EXTRA;S1;")
    assert extra_site_key in result.site_sequences.index
    assert (
        result.site_sequences.at[extra_site_key, "site_sequence"]
        == "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA"
    )
    assert result.dataset_sequences_added == 1
    assert result.dataset_reference_conflict_count == 0
    assert result.conflicts == ()


def test_site_sequence_builder_ignores_matching_dataset_reference_sequences() -> None:
    dataset = _dataset(
        sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
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
            "MAPK14;Y182;": "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
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
        site_sequence_conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    )
    executor = _ExecutorSentinel()

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow(
            interpreter=KinaseWorkflowInterpreter(),
            executor=cast(KinaseWorkflowExecutorContract, executor),
        ).run(request)

    error = exc_info.value
    assert error.seam == "kinase.interpreter.site_sequence_conflict"
    diagnostics = cast(list[dict[str, object]], error.details["conflict_diagnostics"])
    assert len(diagnostics) == 1
    assert diagnostics[0]["site_id"] == "MAPK14;Y182;"
    assert diagnostics[0]["dataset_sequence"] == "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"
    assert diagnostics[0]["reference_sequence"] == "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA"
    assert isinstance(error.next_action, str)
    assert "site_sequence_conflict_policy='prefer_reference'" in error.next_action
    assert "KinaseWorkflowRequest" in error.next_action
    assert executor.called is False


def test_interpreter_default_prefer_reference_records_conflict_diagnostics_in_provenance() -> (
    None
):
    dataset = _dataset(
        sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        }
    )
    result = KinaseWorkflow(interpreter=KinaseWorkflowInterpreter()).run(
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
    assert conflict_rows[0]["selected_sequence"] == "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA"


def test_interpreter_prefer_dataset_selects_dataset_sequence_and_contract_accepts_it() -> (
    None
):
    dataset = _dataset(
        sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAAAAAAAAATAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
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
        site_sequence_conflict_policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    )

    interpreted = KinaseWorkflowInterpreter().run(request)

    mapk14_site_key = site_key_from_display_id("MAPK14;Y182;")
    assert interpreted.site_sequences.at[mapk14_site_key, "site_sequence"] == (
        "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA"
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
            "MAPK14;Y182;": "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        }
    )
    references = _references(
        reference_sequences_by_site={
            "MAPK14;Y182;": "AAAAAAAAAAAAAAAYAAAAAAAAAAAAAAA",
            "GSK3B;S9;": "AAAAAAAAAAAAAAASAAAAAAAAAAAAAAA",
        }
    )
    unnormalised_sequences = references.site_sequences.copy(deep=True)
    unnormalised_sequences.index = pd.Index(
        [" mapk14 ; y182 ", "GSK3B;S9;"],
        name=references.site_sequences.index.name,
    )
    unnormalised_sequences.loc[:, "display_id"] = (
        references.site_sequences.index.astype(str).tolist()
    )

    with pytest.raises(
        WorkflowBoundaryError, match="kinase.contracts.site_sequence_schema"
    ):
        projected_kinase_substrate_map = pd.DataFrame(
            {
                "kinase": ["MAP2K6" for _ in dataset.phospho.index],
                "substrate_site": dataset.phospho.index.astype(str).tolist(),
            }
        )
        ResolvedKinaseWorkflowRequest(
            dataset=dataset,
            references=references,
            kinase_substrate_map=projected_kinase_substrate_map,
            site_sequences=unnormalised_sequences,
            site_identity_map=_site_identity_map(dataset),
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

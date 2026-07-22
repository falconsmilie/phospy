from __future__ import annotations

from typing import NoReturn, cast

import pandas as pd
import pytest

from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflow,
    KinaseWorkflowRequest,
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
)
from phospy.api.configs import (
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
    KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
)
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.workflows.kinase.contracts import ResolvedKinaseWorkflowRequest
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
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

pytestmark = pytest.mark.integration

_DATASET_CONFLICT_SEQUENCE = "AAAAAAASAAAAAAA"
_REFERENCE_CONFLICT_SEQUENCE = "AAAAAAATAAAAAAA"
_MATCHING_SEQUENCE = "AAAAAAAYAAAAAAA"


class _ExecutorMustNotRun:
    called: bool = False

    def run(self, request: ResolvedKinaseWorkflowRequest) -> NoReturn:
        _ = request
        self.called = True
        raise AssertionError("sequence conflict reached kinase executor")


def _dataset(
    *,
    sequences: tuple[str, str] = (_DATASET_CONFLICT_SEQUENCE, _MATCHING_SEQUENCE),
) -> AnalysisReadyPhosphoDataset:
    display_ids = ["GENE1;S10;", "GENE2;Y20;"]
    site_index = site_key_index_from_display_ids(display_ids)
    return trusted_analysis_ready_dataset_from_tables(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0], "sample_b": [1.5, 2.5]},
            index=site_index.copy(),
        ),
        site_metadata=pd.DataFrame(
            {
                "site_key": site_index.astype(str).tolist(),
                "display_id": display_ids,
                **site_key_context_columns(site_index),
                "gene_symbol": ["GENE1", "GENE2"],
                "site": ["S10", "Y20"],
                "protein_id": ["GENE1", "GENE2"],
                "site_sequence": list(sequences),
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _references(
    *,
    sequences: tuple[str, str] = (_REFERENCE_CONFLICT_SEQUENCE, _MATCHING_SEQUENCE),
) -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K1"],
                "substrate_site": ["GENE1;S10;", "GENE2;Y20;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": list(sequences)},
            index=pd.Index(["GENE1;S10;", "GENE2;Y20;"], name="site_id"),
        ),
    )


def _request(
    *,
    policy: object,
    dataset: AnalysisReadyPhosphoDataset | None = None,
    references: ReferenceBundle | None = None,
) -> KinaseWorkflowRequest:
    return KinaseWorkflowRequest(
        dataset=_dataset() if dataset is None else dataset,
        references=_references() if references is None else references,
        scoring_config=KinaseScoringConfig(
            min_substrates=2,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
        prediction_config=KinasePredictionConfig(
            top_k=2,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
            mode="deterministic_ranking",
        ),
        activity_config=None,
        site_sequence_conflict_policy=policy,
    )


def _site_sequence_merge(result: KinaseWorkflowResult) -> dict[str, object]:
    provenance = result.provenance
    assert provenance is not None
    scoring_diagnostics = cast(
        dict[str, object],
        provenance.workflow_parameters["scoring_diagnostics"],
    )
    return cast(dict[str, object], scoring_diagnostics["site_sequence_merge"])


def _single_conflict_row(result: KinaseWorkflowResult) -> dict[str, object]:
    merge = _site_sequence_merge(result)
    rows = cast(list[dict[str, object]], merge["conflict_diagnostics"])
    assert len(rows) == 1
    return rows[0]


@pytest.mark.parametrize(
    "policy",
    (
        KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR,
        KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET,
        KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE,
    ),
)
def test_matching_dataset_reference_sequences_pass_under_all_policies(
    policy: object,
) -> None:
    matching_sequences = (_DATASET_CONFLICT_SEQUENCE, _MATCHING_SEQUENCE)
    dataset = _dataset(sequences=matching_sequences)
    references = _references(sequences=matching_sequences)

    result = KinaseWorkflow().run(
        _request(policy=policy, dataset=dataset, references=references)
    )

    merge = _site_sequence_merge(result)
    assert merge["dataset_reference_conflict_count"] == 0
    assert merge["conflict_diagnostics"] == []


def test_conflicting_sequence_error_policy_fails_before_scoring() -> None:
    executor = _ExecutorMustNotRun()
    site_key = site_key_from_display_id("GENE1;S10;")

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        KinaseWorkflow._with_components(executor=executor).run(
            _request(policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_ERROR)
        )

    assert executor.called is False
    message = str(exc_info.value)
    assert site_key in message
    assert _DATASET_CONFLICT_SEQUENCE in message
    assert _REFERENCE_CONFLICT_SEQUENCE in message
    assert "conflict_policy=error" in message
    diagnostics = cast(
        list[dict[str, object]], exc_info.value.details["conflict_diagnostics"]
    )
    assert diagnostics[0]["site_key"] == site_key
    assert diagnostics[0]["selected_sequence_source"] == "unresolved"


def test_conflicting_sequence_prefer_dataset_records_result_provenance() -> None:
    site_key = site_key_from_display_id("GENE1;S10;")

    result = KinaseWorkflow().run(
        _request(policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET)
    )

    row = _single_conflict_row(result)
    assert row["site_key"] == site_key
    assert row["dataset_sequence"] == _DATASET_CONFLICT_SEQUENCE
    assert row["reference_sequence"] == _REFERENCE_CONFLICT_SEQUENCE
    assert row["selected_sequence"] == _DATASET_CONFLICT_SEQUENCE
    assert row["selected_sequence_source"] == "dataset"
    assert row["policy"] == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_DATASET
    assert "prefer_dataset" in str(row["diagnostic"])
    merge = _site_sequence_merge(result)
    source_rows = cast(list[dict[str, object]], merge["selected_sequence_sources"])
    selected = {str(item["site_key"]): item for item in source_rows}
    assert selected[site_key]["selected_sequence_source"] == "dataset"
    assert selected[site_key]["dataset_sequence"] == _DATASET_CONFLICT_SEQUENCE
    assert selected[site_key]["reference_sequence"] == _REFERENCE_CONFLICT_SEQUENCE


def test_conflicting_sequence_prefer_reference_records_result_provenance() -> None:
    site_key = site_key_from_display_id("GENE1;S10;")

    result = KinaseWorkflow().run(
        _request(policy=KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE)
    )

    row = _single_conflict_row(result)
    assert row["site_key"] == site_key
    assert row["dataset_sequence"] == _DATASET_CONFLICT_SEQUENCE
    assert row["reference_sequence"] == _REFERENCE_CONFLICT_SEQUENCE
    assert row["selected_sequence"] == _REFERENCE_CONFLICT_SEQUENCE
    assert row["selected_sequence_source"] == "reference"
    assert row["policy"] == KINASE_SITE_SEQUENCE_CONFLICT_POLICY_PREFER_REFERENCE
    assert "prefer_reference" in str(row["diagnostic"])
    merge = _site_sequence_merge(result)
    source_rows = cast(list[dict[str, object]], merge["selected_sequence_sources"])
    selected = {str(item["site_key"]): item for item in source_rows}
    assert selected[site_key]["selected_sequence_source"] == "reference"
    assert selected[site_key]["dataset_sequence"] == _DATASET_CONFLICT_SEQUENCE
    assert selected[site_key]["reference_sequence"] == _REFERENCE_CONFLICT_SEQUENCE

from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import (
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from phospy.errors import WorkflowValidationError
from phospy.science.datasets.models import AnalysisReadyPhosphoDataset
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)
from tests.support.unsafe_dataset_states import (
    unsafe_drop_dataset_site_metadata_columns,
    unsafe_set_dataset_site_metadata_columns,
)


def _window(residue: str) -> str:
    return ("A" * 7) + residue + ("A" * 7)


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["GENE1;S10;", "GENE2;T20;"]
    site_index = site_key_index_from_display_ids(display_ids)
    return AnalysisReadyPhosphoDataset(
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
                "site": ["S10", "T20"],
                "protein_id": ["GENE1", "GENE2"],
                "site_sequence": [_window("S"), _window("T")],
            },
            index=site_index.copy(),
        ),
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _kinase_result(
    dataset: AnalysisReadyPhosphoDataset | None = None,
) -> KinaseWorkflowResult:
    resolved_dataset = _dataset() if dataset is None else dataset
    site_index = resolved_dataset.phospho.index.copy()
    score_matrix = pd.DataFrame({"K1": [1.0, 2.0]}, index=site_index.copy())
    prediction_matrix = pd.DataFrame({"K1": [0.8, 0.9]}, index=site_index.copy())
    return KinaseWorkflowResult(
        dataset=resolved_dataset,
        references=ReferenceBundle(
            organism=Organism.RAT,
            kinase_substrate_map=pd.DataFrame(
                {
                    "kinase": ["K1", "K1"],
                    "substrate_site": ["GENE1;S10;", "GENE2;T20;"],
                }
            ),
            site_sequences=pd.DataFrame(
                {"site_sequence": [_window("S"), _window("T")]},
                index=pd.Index(["GENE1;S10;", "GENE2;T20;"], name="site_id"),
            ),
        ),
        scoring_result=KinaseScoringResult._from_owned(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix,
        ),
        prediction_result=KinasePredictionResult._from_owned(
            pred_mat=prediction_matrix
        ),
        activity_result=None,
    )


def _request(
    dataset: AnalysisReadyPhosphoDataset | None = None,
) -> SignalomeWorkflowRequest:
    return SignalomeWorkflowRequest(
        kinase_result=_kinase_result(dataset),
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            reference_context_compatibility_policy=(
                ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
            ),
        ),
    )


def test_signalome_accepts_valid_centered_sequence_context() -> None:
    request = _request()

    validated = SignalomeWorkflowValidator().run(request)

    assert validated is request


def test_signalome_requires_site_sequence_column() -> None:
    dataset = _dataset()
    unsafe_drop_dataset_site_metadata_columns(dataset, "site_sequence")

    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(_request(dataset))

    message = str(exc_info.value)
    assert "signalome workflow request" in message
    assert "site_sequence" in message


@pytest.mark.parametrize(
    ("sequence", "expected"),
    [
        ("", "missing or blank"),
        ("AAAAAAAA", "must be odd length"),
        ("AAAAAAATAAAAAAA", "centre residue must match the site token residue"),
        ("AAAAAAA*AAAAAAA", "unsupported_characters='*'"),
    ],
)
def test_signalome_rejects_invalid_centered_sequence_context(
    sequence: str,
    expected: str,
) -> None:
    dataset = _dataset()
    unsafe_set_dataset_site_metadata_columns(
        dataset,
        {"site_sequence": [sequence, _window("T")]},
    )

    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(_request(dataset))

    message = str(exc_info.value)
    assert "requires centred sequence context" in message
    assert str(dataset.site_metadata.index[0]) in message
    assert expected in message

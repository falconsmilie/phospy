from __future__ import annotations

import pandas as pd
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    DatasetBuildRequest,
    Organism,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
)
from phospy.errors import PhosPyInputError, WorkflowValidationError
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.science.datasets.builders.interpreter import DatasetBuildRequestInterpreter
from phospy.science.references.models import ReferenceBundle
from phospy.workflows.signalome.constants import (
    SIGNALOME_INTERPRETER_SITE_ALIGNMENT_SEAM,
)
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter
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


def test_workflow_boundary_error_supports_message_only_construction() -> None:
    error = WorkflowBoundaryError("workflow boundary message")

    assert str(error) == "workflow boundary message"
    assert error.seam is None
    assert error.next_action is None
    assert error.details == {}


def test_workflow_boundary_error_exposes_structured_diagnostics() -> None:
    error = WorkflowBoundaryError(
        seam="kinase.interpreter.reference_coverage",
        next_action="use compatible references",
        details={"dataset_sites": 2, "overlap_sites": 0},
        message_prefix="kinase workflow boundary validation failed",
    )

    assert error.seam == "kinase.interpreter.reference_coverage"
    assert error.next_action == "use compatible references"
    assert error.details == {"dataset_sites": 2, "overlap_sites": 0}
    assert "seam=kinase.interpreter.reference_coverage" in str(error)
    assert "dataset_sites=2" in str(error)
    assert "next_action=use compatible references" in str(error)


def test_workflow_boundary_error_copies_details_mapping() -> None:
    details = {"shared_sites": 0}
    error = WorkflowBoundaryError(
        seam=SIGNALOME_INTERPRETER_SITE_ALIGNMENT_SEAM,
        next_action="align score and prediction matrices",
        details=details,
    )

    details["shared_sites"] = 99

    assert error.details == {"shared_sites": 0}


def _dataset(*, with_protein_id: bool) -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "AKT1;T308;"]
    index = _site_index()
    site_metadata = pd.DataFrame(
        {
            "site_key": index.tolist(),
            "display_id": display_ids,
            **site_key_context_columns(index),
            "gene_symbol": ["MAPK14", "AKT1"],
            "site": ["Y182", "T308"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "T308"]
            ],
        },
        index=index,
    )
    if with_protein_id:
        site_metadata.loc[:, "protein_id"] = ["P53778", "P31749"]
    return AnalysisReadyPhosphoDataset(
        phospho=pd.DataFrame(
            {"sample_a": [1.0, 2.0], "sample_b": [3.0, 4.0]},
            index=index,
        ),
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _site_index() -> pd.Index:
    return site_key_index_from_display_ids(
        ["MAPK14;Y182;", "AKT1;T308;"],
        protein_namespace="gene_symbol",
    )


def _bundle() -> ReferenceBundle:
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "AKT1"],
                "substrate_site": ["MAPK14;Y182;", "AKT1;T308;"],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "B" * 31]},
            index=pd.Index(["MAPK14;Y182;", "AKT1;T308;"], name="site_id"),
        ),
    )


def _kinase_result(
    *,
    prediction_matrix: pd.DataFrame,
    score_matrix: pd.DataFrame,
    dataset: AnalysisReadyPhosphoDataset | None = None,
) -> KinaseWorkflowResult:
    return KinaseWorkflowResult(
        dataset=_dataset(with_protein_id=True) if dataset is None else dataset,
        references=_bundle(),
        scoring_result=KinaseScoringResult(profile_scores=score_matrix),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
    )


def _request(
    *,
    prediction_matrix: pd.DataFrame,
    score_matrix: pd.DataFrame,
    dataset: AnalysisReadyPhosphoDataset | None = None,
) -> SignalomeWorkflowRequest:
    return SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
            dataset=dataset,
        ),
        config=build_signalome_config(),
    )


def test_signalome_boundary_wraps_non_numeric_downstream_score_matrix() -> None:
    request = _request(
        prediction_matrix=pd.DataFrame(
            {"MAP2K6": [0.8, 0.7]},
            index=_site_index(),
        ),
        score_matrix=pd.DataFrame(
            {"MAP2K6": [1.0, 2.0]},
            index=_site_index(),
        ),
    )
    bad_score_matrix = pd.DataFrame(
        {"MAP2K6": ["abc", "2.0"]},
        index=_site_index(),
    )
    object.__setattr__(
        request.kinase_result.scoring_result, "_profile_scores", bad_score_matrix
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowInterpreter().run(request)

    message = str(exc_info.value)
    assert (
        "signalome.score_matrix_conversion failed while converting downstream score matrix to float"
        in message
    )
    assert "Original error: ValueError:" in message
    assert (
        "Next action: ensure kinase scoring outputs contain numeric finite values"
        in message
    )


def test_signalome_boundary_wraps_non_numeric_prediction_matrix() -> None:
    request = _request(
        prediction_matrix=pd.DataFrame(
            {"MAP2K6": [0.8, 0.7]},
            index=_site_index(),
        ),
        score_matrix=pd.DataFrame(
            {"MAP2K6": [1.0, 2.0]},
            index=_site_index(),
        ),
    )
    bad_prediction_matrix = pd.DataFrame(
        {"MAP2K6": ["bad", "0.7"]},
        index=_site_index(),
    )
    object.__setattr__(
        request.kinase_result.prediction_result, "_pred_mat", bad_prediction_matrix
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowInterpreter().run(request)

    message = str(exc_info.value)
    assert (
        "signalome.prediction_matrix_conversion failed while converting prediction matrix to float"
        in message
    )
    assert "Original error: ValueError:" in message
    assert (
        "Next action: ensure kinase prediction outputs contain numeric finite values"
        in message
    )


def test_signalome_validator_rejects_missing_required_site_metadata_column() -> None:
    dataset = _dataset(with_protein_id=False)
    request = _request(
        prediction_matrix=pd.DataFrame(
            {"MAP2K6": [0.8, 0.7]},
            index=_site_index(),
        ),
        score_matrix=pd.DataFrame(
            {"MAP2K6": [1.0, 2.0]},
            index=_site_index(),
        ),
        dataset=dataset,
    )

    with pytest.raises(WorkflowValidationError) as exc_info:
        SignalomeWorkflowValidator().run(request)

    message = str(exc_info.value)
    assert "is missing required columns: protein_id" in message
    assert (
        "Signalome execution requires an explicit site_metadata.protein_id column"
        in message
    )


def test_signalome_boundary_rejects_misaligned_prediction_and_scoring_indices() -> None:
    request = _request(
        prediction_matrix=pd.DataFrame(
            {"MAP2K6": [0.8, 0.7]},
            index=site_key_index_from_display_ids(["A;S1;", "B;S2;"]),
        ),
        score_matrix=pd.DataFrame(
            {"MAP2K6": [1.0, 2.0]},
            index=site_key_index_from_display_ids(["C;S3;", "D;S4;"]),
        ),
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowInterpreter().run(request)

    error = exc_info.value
    message = str(error)
    assert error.seam == SIGNALOME_INTERPRETER_SITE_ALIGNMENT_SEAM
    assert "shared_sites=0" in message
    assert "next_action=" in message


def test_signalome_boundary_rejects_infinite_downstream_scores() -> None:
    request = _request(
        prediction_matrix=pd.DataFrame(
            {"MAP2K6": [0.8, 0.7]},
            index=_site_index(),
        ),
        score_matrix=pd.DataFrame(
            {"MAP2K6": [1.0, 2.0]},
            index=_site_index(),
        ),
    )
    bad_score_matrix = pd.DataFrame(
        {"MAP2K6": [float("inf"), 2.0]},
        index=_site_index(),
    )
    object.__setattr__(
        request.kinase_result.scoring_result, "_profile_scores", bad_score_matrix
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowInterpreter().run(request)

    message = str(exc_info.value)
    assert "signalome workflow boundary validation failed" in message
    assert "infinite_score_entries=1" in message
    assert "next_action=" in message


def test_builder_interpreter_wraps_invalid_matrix_value_failures() -> None:
    class _FailingNormalizer:
        def run(self, **_: object) -> object:
            raise ValueError("could not convert string to float: 'abc'")

    request = DatasetBuildRequest(
        phospho=pd.DataFrame(
            {"sample_a": [1.0]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
        site_metadata=pd.DataFrame(
            {"gene_symbol": ["MAPK14"], "site": ["Y182"], "site_sequence": ["A" * 31]},
            index=pd.Index(["MAPK14;Y182;"], name="site_id"),
        ),
    )

    with pytest.raises(PhosPyInputError) as exc_info:
        DatasetBuildRequestInterpreter(normalizer=_FailingNormalizer()).run(request)

    message = str(exc_info.value)
    assert (
        "dataset_builder.normalization failed while normalizing input indices and metadata column conventions"
        in message
    )
    assert (
        "Original error: ValueError: could not convert string to float: 'abc'"
        in message
    )
    assert (
        "Next action: ensure phospho/site_metadata/sample_metadata/total tables use supported rectangular DataFrame shapes"
        in message
    )

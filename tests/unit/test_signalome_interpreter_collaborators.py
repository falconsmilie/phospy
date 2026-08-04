from __future__ import annotations

import pandas as pd
import pandas.testing as pdt
import pytest

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    SignalomeWorkflowRequest,
)
from phospy.api.configs import (
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
)
from phospy.api.results import KinasePredictionResult, KinaseScoringResult
from phospy.errors.workflows import WorkflowBoundaryError
from phospy.workflows.signalome.alignment_diagnostics import (
    SignalomeAlignmentDiagnosticsBuilder,
)
from phospy.workflows.signalome.constants import (
    SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM,
)
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter
from phospy.workflows.signalome.matrix_alignment import SignalomeMatrixAligner
from phospy.workflows.signalome.protein_resolution import SignalomeProteinResolver
from phospy.workflows.signalome.score_preconditioning import (
    SignalomeScorePreconditioner,
)
from tests.support.analysis_ready_dataset_factories import (
    trusted_analysis_ready_dataset_from_tables,
)
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)


def _dataset(
    *,
    site_ids: list[str],
    protein_ids: list[str] | None = None,
) -> AnalysisReadyPhosphoDataset:
    if protein_ids is None:
        protein_ids = [site_id.split(";", 1)[0] for site_id in site_ids]
    site_index = site_key_index_from_display_ids(site_ids)
    sites = [site_id.split(";")[1] for site_id in site_ids]
    phospho = pd.DataFrame(
        {
            "sample_a": [float(i + 1) for i in range(len(site_ids))],
            "sample_b": [float(i + 2) for i in range(len(site_ids))],
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": site_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": [site_id.split(";", 1)[0] for site_id in site_ids],
            "site": sites,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "protein_id": protein_ids,
        },
        index=site_index.copy(),
    )
    return trusted_analysis_ready_dataset_from_tables(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _bundle(site_ids: list[str]) -> ReferenceBundle:
    unique_sites = pd.Index(site_ids, name="site_id")
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["K1"], "substrate_site": [str(unique_sites[0])]}
        ),
        site_sequences=pd.DataFrame(
            {
                "site_sequence": [
                    ("A" * 15)
                    + str(site_id).split(";")[1].strip().upper()[0]
                    + ("A" * 15)
                    for site_id in unique_sites
                ]
            },
            index=unique_sites,
        ),
    )


def _matrix(
    *,
    values: list[list[object]],
    site_ids: list[str],
    kinases: list[str],
) -> pd.DataFrame:
    return pd.DataFrame(
        values,
        index=_site_index(site_ids),
        columns=pd.Index(kinases, name="kinase"),
    )


def _site_index(site_ids: list[str]) -> pd.Index:
    if all(str(site_id).startswith("phospy:v1|") for site_id in site_ids):
        return pd.Index(site_ids, name="site_key")
    if all(str(site_id).count(";") >= 2 for site_id in site_ids):
        return site_key_index_from_display_ids(site_ids)
    return pd.Index(site_ids, name="site_id")


def _site_keys(site_ids: list[str]) -> list[str]:
    return site_key_index_from_display_ids(site_ids).astype(str).tolist()


def _request(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    prediction_matrix: pd.DataFrame,
    profile_scores: pd.DataFrame,
    rank_weighted_scores: pd.DataFrame | None = None,
) -> SignalomeWorkflowRequest:
    display_site_ids = dataset.site_metadata.loc[:, "display_id"].astype(str).tolist()
    return SignalomeWorkflowRequest(
        kinase_result=KinaseWorkflowResult(
            dataset=dataset,
            references=_bundle(site_ids=display_site_ids),
            scoring_result=KinaseScoringResult(
                profile_scores=profile_scores,
                rank_weighted_fusion_scores=rank_weighted_scores,
            ),
            prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
            activity_result=None,
        ),
        config=build_signalome_config(),
    )


def test_matrix_aligner_preserves_shared_site_and_kinase_ordering() -> None:
    aligner = SignalomeMatrixAligner()
    aligned = aligner.run(
        dataset_sites=pd.Index(["S3", "S1", "S2"]),
        prediction_matrix=_matrix(
            values=[[1.0, 2.0, 3.0], [1.1, 2.1, 3.1], [1.2, 2.2, 3.2], [9.0, 9.0, 9.0]],
            site_ids=["S2", "S1", "S3", "S4"],
            kinases=["K2", "K1", "K4"],
        ),
        downstream_score_matrix=_matrix(
            values=[[4.0, 5.0, 6.0], [4.1, 5.1, 6.1], [4.2, 5.2, 6.2], [7.0, 8.0, 9.0]],
            site_ids=["S1", "S3", "S2", "S5"],
            kinases=["K1", "K2", "K3"],
        ),
        downstream_score_source="rank_weighted_fusion_scores",
    )

    assert aligned.aligned_site_index.tolist() == ["S3", "S1", "S2"]
    assert aligned.aligned_kinase_index.tolist() == ["K2", "K1"]
    assert aligned.aligned_prediction_matrix.index.tolist() == ["S3", "S1", "S2"]
    assert aligned.aligned_prediction_matrix.columns.tolist() == ["K2", "K1"]
    assert aligned.aligned_downstream_score_matrix.index.tolist() == ["S3", "S1", "S2"]
    assert aligned.aligned_downstream_score_matrix.columns.tolist() == ["K2", "K1"]


def test_matrix_aligner_rejects_duplicate_labels_with_wrapped_boundary_error() -> None:
    aligner = SignalomeMatrixAligner()
    with pytest.raises(WorkflowBoundaryError) as exc_info:
        aligner.run(
            dataset_sites=pd.Index(["S1", "S2"]),
            prediction_matrix=_matrix(
                values=[[0.9], [0.8]],
                site_ids=["S1", "S1"],
                kinases=["K1"],
            ),
            downstream_score_matrix=_matrix(
                values=[[1.0], [2.0]],
                site_ids=["S1", "S2"],
                kinases=["K1"],
            ),
            downstream_score_source="rank_weighted_fusion_scores",
        )

    error = exc_info.value
    message = str(error)
    assert error.seam == "signalome.interpreter.prediction_matrix_conversion"
    assert "converting prediction matrix to float" in message
    assert error.details["field_name"].endswith("prediction_result.pred_mat")
    assert error.details["original_error_type"] == "WorkflowValidationError"


def test_matrix_aligner_wraps_numeric_coercion_failures_as_actionable_boundary_errors() -> (
    None
):
    aligner = SignalomeMatrixAligner()
    with pytest.raises(WorkflowBoundaryError) as exc_info:
        aligner.run(
            dataset_sites=pd.Index(["S1", "S2"]),
            prediction_matrix=_matrix(
                values=[["bad"], [0.8]],
                site_ids=["S1", "S2"],
                kinases=["K1"],
            ),
            downstream_score_matrix=_matrix(
                values=[[1.0], [2.0]],
                site_ids=["S1", "S2"],
                kinases=["K1"],
            ),
            downstream_score_source="rank_weighted_fusion_scores",
        )

    message = str(exc_info.value)
    assert (
        "signalome.prediction_matrix_conversion failed while converting prediction matrix to float"
        in message
    )
    assert (
        "Next action: ensure kinase prediction outputs contain numeric finite values"
        in message
    )


def test_score_preconditioner_strict_policy_errors_on_all_missing_rows() -> None:
    preconditioner = SignalomeScorePreconditioner()
    with pytest.raises(WorkflowBoundaryError) as exc_info:
        preconditioner.run(
            score_matrix=_matrix(
                values=[[1.0, 2.0], [float("nan"), float("nan")]],
                site_ids=["S1", "S2"],
                kinases=["K1", "K2"],
            ),
            policy=SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
        )

    error = exc_info.value
    assert error.seam == SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM
    assert error.details["dropped_all_missing_row_count"] == 1
    assert error.details["retained_row_count"] == 1


def test_score_preconditioner_permissive_policy_drops_and_reports_rows() -> None:
    preconditioner = SignalomeScorePreconditioner()
    result = preconditioner.run(
        score_matrix=_matrix(
            values=[[1.0, 2.0], [float("nan"), float("nan")], [3.0, 4.0]],
            site_ids=["S1", "S2", "S3"],
            kinases=["K1", "K2"],
        ),
        policy=SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    )

    assert result.downstream_score_matrix.index.tolist() == ["S1", "S3"]
    assert result.diagnostics.input_row_count == 3
    assert result.diagnostics.dropped_all_missing_row_count == 1
    assert result.diagnostics.retained_row_count == 2
    assert (
        result.diagnostics.policy
        == SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )


def test_score_preconditioner_rejects_infinite_values() -> None:
    preconditioner = SignalomeScorePreconditioner()
    with pytest.raises(WorkflowBoundaryError) as exc_info:
        preconditioner.run(
            score_matrix=_matrix(
                values=[[float("inf"), 2.0]],
                site_ids=["S1"],
                kinases=["K1", "K2"],
            ),
            policy=SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
        )

    error = exc_info.value
    assert error.seam == SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM
    assert error.details["infinite_score_entries"] == 1


def test_protein_resolver_maps_retained_sites_to_explicit_protein_group_ids() -> None:
    resolver = SignalomeProteinResolver()
    dataset = _dataset(
        site_ids=["P1;S1;", "P2;S2;", "P3;S3;"],
        protein_ids=["PROT1", "PROT2", "PROT3"],
    )

    resolved = resolver.run(
        site_metadata=dataset.site_metadata,
        site_index=site_key_index_from_display_ids(["P3;S3;", "P1;S1;"]),
        removed_by_score_preconditioning_count=1,
    )

    assert resolved.index.tolist() == _site_keys(["P3;S3;", "P1;S1;"])
    assert resolved.tolist() == ["PROT3", "PROT1"]
    assert resolved.name == "protein_group_id"


def test_alignment_diagnostics_builder_reports_expected_counts_and_reasons() -> None:
    diagnostics = SignalomeAlignmentDiagnosticsBuilder().run(
        dataset_sites=pd.Index(["S1", "S2", "S3"]),
        prediction_sites=pd.Index(["S2", "S3", "S4"]),
        score_sites=pd.Index(["S3", "S2", "S5"]),
        shared_sites=pd.Index(["S2", "S3"]),
        retained_sites=pd.Index(["S3"]),
        prediction_kinases=pd.Index(["K2", "K1", "K3"]),
        score_kinases=pd.Index(["K1", "K4", "K2"]),
        shared_kinases=pd.Index(["K2", "K1"]),
        interpreted_protein_sites=pd.Index(["S2", "S3"]),
        retained_protein_sites=pd.Index(["S3"]),
    )

    assert diagnostics.dataset_sites.dropped_reasons == {
        "missing_from_prediction_scores": 1,
        "missing_from_downstream_scores": 0,
        "removed_by_score_preconditioning": 1,
        "removed_by_validation_policy": 0,
    }
    assert diagnostics.prediction_score_sites.dropped_reasons == {
        "missing_from_dataset": 1,
        "missing_from_downstream_scores": 0,
        "removed_by_score_preconditioning": 1,
        "removed_by_validation_policy": 0,
    }
    assert diagnostics.downstream_score_sites.dropped_reasons == {
        "missing_from_dataset": 1,
        "missing_from_prediction_scores": 0,
        "removed_by_score_preconditioning": 1,
        "removed_by_validation_policy": 0,
    }
    assert diagnostics.kinases.dropped_reasons == {
        "missing_from_prediction_scores": 1,
        "missing_from_downstream_scores": 1,
        "missing_kinase_support": 0,
    }
    assert diagnostics.protein_group_ids.dropped_reasons == {
        "removed_by_score_preconditioning": 1,
        "missing_protein_group_id": 0,
        "removed_by_validation_policy": 0,
    }


def test_interpreter_run_preserves_expected_interpreted_output_for_valid_fixture() -> (
    None
):
    dataset = _dataset(
        site_ids=["P1;S1;", "P2;S2;", "P3;S3;"],
        protein_ids=["PROT1", "PROT2", "PROT3"],
    )
    prediction_matrix = _matrix(
        values=[[0.9, 0.1], [0.2, 0.8], [0.3, 0.7]],
        site_ids=["P1;S1;", "P2;S2;", "P3;S3;"],
        kinases=["K1", "K2"],
    ).astype(float)
    profile_scores = _matrix(
        values=[[5.0, 5.0], [5.0, 5.0], [5.0, 5.0]],
        site_ids=["P1;S1;", "P2;S2;", "P3;S3;"],
        kinases=["K1", "K2"],
    ).astype(float)
    rank_weighted_scores = _matrix(
        values=[[1.0, 2.0], [float("nan"), float("nan")], [3.0, 4.0]],
        site_ids=["P1;S1;", "P2;S2;", "P3;S3;"],
        kinases=["K1", "K2"],
    ).astype(float)
    request = _request(
        dataset=dataset,
        prediction_matrix=prediction_matrix,
        profile_scores=profile_scores,
        rank_weighted_scores=rank_weighted_scores,
    )
    request = SignalomeWorkflowRequest(
        kinase_result=request.kinase_result,
        config=build_signalome_config(score_preconditioning_policy="allow_and_report"),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)

    assert interpreted.downstream_score_source == "rank_weighted_fusion_scores"
    retained_site_keys = _site_keys(["P1;S1;", "P3;S3;"])
    assert interpreted.downstream_score_matrix.index.tolist() == retained_site_keys
    assert interpreted.prediction_matrix.index.tolist() == retained_site_keys
    assert interpreted.site_to_protein_group_id.index.tolist() == retained_site_keys
    assert interpreted.site_to_protein_group_id.tolist() == ["PROT1", "PROT3"]
    pdt.assert_index_equal(
        interpreted.downstream_score_matrix.columns,
        interpreted.prediction_matrix.columns,
    )
    assert interpreted.score_preconditioning_diagnostics.input_row_count == 3
    assert (
        interpreted.score_preconditioning_diagnostics.dropped_all_missing_row_count == 1
    )
    assert interpreted.score_preconditioning_diagnostics.retained_row_count == 2

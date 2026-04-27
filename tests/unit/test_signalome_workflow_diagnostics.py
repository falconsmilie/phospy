from __future__ import annotations

import warnings

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyPhosphoDataset,
)
from phospy.api import (
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    SignalomeConfig,
    SignalomeWorkflowRequest,
)
from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
    SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
    SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
    SIGNALOME_CLUSTERING_BACKEND_APPROXIMATE,
    SIGNALOME_CLUSTERING_BACKEND_EXACT,
    SIGNALOME_MAX_EXACT_CLUSTER_TREE_SITES_DEFAULT,
    SIGNALOME_MAX_FULL_CORRELATION_SITES_DEFAULT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    SignalomeWorkflowResult,
)
from phospy.errors import (
    SignalomeScaleError,
    WorkflowBoundaryError,
    WorkflowValidationError,
)
from phospy.errors.workflows import WorkflowStageError
from phospy.signalomes.clustering import (
    MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD,
    SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY,
)
from phospy.signalomes.constants import (
    EXPANDED_SIGNALOME_ROW_KIND_SUMMARY,
    SITE_CLUSTER_COLUMN,
)
from phospy.signalomes.models import (
    KinaseNetwork,
    SignalomeAssignments,
    SignalomeModules,
)
from phospy.workflows.signalome.constants import (
    SIGNALOME_EXECUTOR_CONTEXT_TABLES_SEAM,
    SIGNALOME_EXECUTOR_EXPANDED_SIGNALOME_SEAM,
    SIGNALOME_EXECUTOR_KINASE_SUPPORT_SEAM,
    SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM,
    SIGNALOME_EXECUTOR_NETWORK_SEAM,
    SIGNALOME_INTERPRETER_KINASE_OVERLAP_SEAM,
    SIGNALOME_INTERPRETER_PROTEIN_MAPPING_SEAM,
    SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM,
    SIGNALOME_INTERPRETER_SITE_ALIGNMENT_SEAM,
    SIGNALOME_PROTEIN_RESOLUTION_SOURCE_SITE_METADATA,
)
from phospy.workflows.signalome.contracts import (
    ResolvedSignalomeExecutionConfig,
    ResolvedSignalomeWorkflowRequest,
)
from phospy.workflows.signalome.executor import SignalomeWorkflowExecutor
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)


def _dataset(
    *,
    site_ids: list[str],
    gene_symbols: list[str] | None = None,
    protein_ids: list[str] | None = None,
) -> AnalysisReadyPhosphoDataset:
    if gene_symbols is None:
        gene_symbols = [str(site_id).split(";", 1)[0] for site_id in site_ids]
    if protein_ids is None:
        protein_ids = [str(site_id).split(";", 1)[0].strip() for site_id in site_ids]
    phospho = pd.DataFrame(
        {
            "sample_a": [float(index + 1) for index in range(len(site_ids))],
            "sample_b": [float(index + 2) for index in range(len(site_ids))],
        },
        index=site_ids,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": gene_symbols,
            "site": [f"S{index + 1}" for index in range(len(site_ids))],
            "site_sequence": ["A" * 31 for _ in site_ids],
            "protein_id": protein_ids,
        },
        index=site_ids,
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        intensity_scale_state=supported_linear_intensity_scale_state(
            has_total_matrix=False
        ),
        processing_state=supported_linear_processing_state(has_total_matrix=False),
    )


def _bundle(site_ids: list[str]) -> ReferenceBundle:
    unique_sites = pd.Index([str(site_id) for site_id in site_ids], name="site_id")
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {"kinase": ["MAP2K6"], "substrate_site": [str(unique_sites[0])]}
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31 for _ in unique_sites]},
            index=unique_sites,
        ),
    )


def _matrix(
    *,
    values: list[list[float]],
    site_ids: list[str],
    kinases: list[str],
) -> pd.DataFrame:
    return pd.DataFrame(
        values,
        index=pd.Index(site_ids, name="site_id"),
        columns=pd.Index(kinases, name="kinase"),
        dtype=float,
    )


def _kinase_result(
    *,
    dataset: AnalysisReadyPhosphoDataset,
    prediction_matrix: pd.DataFrame,
    score_matrix: pd.DataFrame,
    combined_score_matrix: pd.DataFrame | None = None,
) -> KinaseWorkflowResult:
    return KinaseWorkflowResult(
        dataset=dataset,
        references=_bundle(site_ids=dataset.phospho.index.astype(str).tolist()),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=(
                score_matrix if combined_score_matrix is None else combined_score_matrix
            ),
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


def _execution_config(config: SignalomeConfig) -> ResolvedSignalomeExecutionConfig:
    return ResolvedSignalomeExecutionConfig(
        substrate_support_cutoff=float(config.substrate_support_cutoff),
        network_correlation_threshold=float(config.network_correlation_threshold),
        network_policy=config.network_policy,
        assignment_policy=config.assignment_policy,
        score_preconditioning_policy=config.score_preconditioning_policy,
        module_selection_primary_threshold=float(
            config.module_selection_primary_correlation_threshold
        ),
        module_selection_fallback_threshold=float(
            config.module_selection_fallback_correlation_threshold
        ),
        module_selection_max_clusters=int(config.module_selection_max_clusters),
        cluster_tree_backend=config.cluster_tree_backend,
        candidate_scoring_backend=config.candidate_scoring_backend,
        max_exact_cluster_tree_sites=int(config.max_exact_cluster_tree_sites),
        max_full_correlation_sites=int(config.max_full_correlation_sites),
        clustering_backend_alias=config.clustering_backend,
        requested_module_count=(
            None if config.module_count is None else int(config.module_count)
        ),
    )


def _interpreted_request_for_context_failure() -> ResolvedSignalomeWorkflowRequest:
    dataset = _dataset(site_ids=["P1;S1;", "P2;S2;", "P3;S3;"])
    prediction_matrix = _matrix(
        values=[
            [0.95, 0.1],
            [0.1, 0.95],
            [0.8, 0.7],
        ],
        site_ids=["P1;S1;", "P2;S2;", "P3;S3;"],
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, 0.0],
            [0.0, 1.0],
            [0.8, 0.7],
        ],
        site_ids=["P1;S1;", "P2;S2;", "P3;S3;"],
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )
    return SignalomeWorkflowInterpreter().run(request)


def test_boundary_error_reports_no_usable_site_alignment_counts() -> None:
    dataset = _dataset(site_ids=["P1;S1;", "P2;S2;"])
    prediction_matrix = _matrix(
        values=[[0.9], [0.8]],
        site_ids=["X1;S1;", "X2;S2;"],
        kinases=["K1"],
    )
    score_matrix = _matrix(
        values=[[1.0], [2.0]],
        site_ids=["X1;S1;", "X2;S2;"],
        kinases=["K1"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowInterpreter().run(request)

    error = exc_info.value
    message = str(error)
    assert error.seam == SIGNALOME_INTERPRETER_SITE_ALIGNMENT_SEAM
    assert error.next_action is not None
    assert error.details["dataset_sites"] == 2
    assert error.details["prediction_sites"] == 2
    assert error.details["score_sites"] == 2
    assert error.details["shared_sites"] == 0
    assert f"seam={SIGNALOME_INTERPRETER_SITE_ALIGNMENT_SEAM}" in message
    assert "dataset_sites=2" in message
    assert "prediction_sites=2" in message
    assert "score_sites=2" in message
    assert "shared_sites=0" in message
    assert "next_action=" in message


def test_boundary_error_reports_no_overlapping_kinase_set_counts() -> None:
    dataset = _dataset(site_ids=["P1;S1;", "P2;S2;"])
    prediction_matrix = _matrix(
        values=[[0.9, 0.1], [0.2, 0.8]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[[1.0, 2.0], [3.0, 4.0]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["A1", "A2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowInterpreter().run(request)

    message = str(exc_info.value)
    assert f"seam={SIGNALOME_INTERPRETER_KINASE_OVERLAP_SEAM}" in message
    assert "prediction_kinases=2" in message
    assert "score_kinases=2" in message
    assert "shared_kinases=0" in message
    assert "next_action=" in message


def test_boundary_error_reports_unusable_protein_mapping_counts() -> None:
    base_dataset = _dataset(
        site_ids=["MAPK14;S1;", "GSK3B;S2;"],
        gene_symbols=["MAPK14", "GSK3B"],
        protein_ids=["P28482", "Q9Y243"],
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=base_dataset.phospho,
        site_metadata=base_dataset.site_metadata.drop(columns=["protein_id"]),
        sample_metadata=base_dataset.sample_metadata,
        total=base_dataset.total,
        organism=base_dataset.organism,
        intensity_scale_state=base_dataset.intensity_scale_state,
        processing_state=base_dataset.processing_state,
    )
    prediction_matrix = _matrix(
        values=[[0.9], [0.8]],
        site_ids=["MAPK14;S1;", "GSK3B;S2;"],
        kinases=["K1"],
    )
    score_matrix = _matrix(
        values=[[1.0], [2.0]],
        site_ids=["MAPK14;S1;", "GSK3B;S2;"],
        kinases=["K1"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowInterpreter().run(request)

    message = str(exc_info.value)
    assert f"seam={SIGNALOME_INTERPRETER_PROTEIN_MAPPING_SEAM}" in message
    assert (
        f"protein_resolution_source={SIGNALOME_PROTEIN_RESOLUTION_SOURCE_SITE_METADATA}"
    ) in message
    assert "interpreted_sites=2" in message
    assert "resolved_protein_sites=0" in message
    assert "unresolved_protein_sites=2" in message
    assert "next_action=" in message


def test_interpreter_uses_explicit_site_metadata_protein_id_when_present() -> None:
    dataset = _dataset(
        site_ids=["MAPK14;S1;", "MAPK14;S2;"],
        gene_symbols=["MAPK14", "MAPK14"],
        protein_ids=["P28482-1", "P28482-2"],
    )
    prediction_matrix = _matrix(
        values=[[0.9], [0.8]],
        site_ids=["MAPK14;S1;", "MAPK14;S2;"],
        kinases=["K1"],
    )
    score_matrix = _matrix(
        values=[[1.0], [2.0]],
        site_ids=["MAPK14;S1;", "MAPK14;S2;"],
        kinases=["K1"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)
    assert interpreted.site_to_protein.tolist() == ["P28482-1", "P28482-2"]


def test_interpreter_does_not_fallback_to_site_id_prefix_when_protein_id_column_missing() -> (
    None
):
    base_dataset = _dataset(
        site_ids=["MAPK14;S1;", "MAPK14;S2;"],
        gene_symbols=["MAPK14", "MAPK14"],
        protein_ids=["P28482-1", "P28482-2"],
    )
    dataset = AnalysisReadyPhosphoDataset(
        phospho=base_dataset.phospho,
        site_metadata=base_dataset.site_metadata.drop(columns=["protein_id"]),
        sample_metadata=base_dataset.sample_metadata,
        total=base_dataset.total,
        organism=base_dataset.organism,
        intensity_scale_state=base_dataset.intensity_scale_state,
        processing_state=base_dataset.processing_state,
    )
    prediction_matrix = _matrix(
        values=[[0.9], [0.8]],
        site_ids=["MAPK14;S1;", "MAPK14;S2;"],
        kinases=["K1"],
    )
    score_matrix = _matrix(
        values=[[1.0], [2.0]],
        site_ids=["MAPK14;S1;", "MAPK14;S2;"],
        kinases=["K1"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowInterpreter().run(request)

    message = str(exc_info.value)
    assert f"seam={SIGNALOME_INTERPRETER_PROTEIN_MAPPING_SEAM}" in message
    assert (
        f"protein_resolution_source={SIGNALOME_PROTEIN_RESOLUTION_SOURCE_SITE_METADATA}"
    ) in message
    assert "resolved_protein_sites=0" in message


def test_interpreter_prefers_rank_weighted_fusion_scores_for_downstream_signalome_matrix() -> (
    None
):
    dataset = _dataset(site_ids=["P1;S1;", "P2;S2;"])
    prediction_matrix = _matrix(
        values=[[0.9, 0.2], [0.1, 0.8]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1", "K2"],
    )
    profile_scores = _matrix(
        values=[[0.1, 0.9], [0.8, 0.2]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1", "K2"],
    )
    rank_weighted_fusion_scores = _matrix(
        values=[[0.7, 0.4], [0.3, 0.6]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=profile_scores,
            combined_score_matrix=rank_weighted_fusion_scores,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)
    pd.testing.assert_frame_equal(
        interpreted.downstream_score_matrix,
        rank_weighted_fusion_scores,
        check_dtype=False,
    )
    assert interpreted.downstream_score_source == "rank_weighted_fusion_scores"


def test_interpreter_resolves_execution_config_defaults_for_executor() -> None:
    dataset = _dataset(site_ids=["P1;S1;", "P2;S2;"])
    prediction_matrix = _matrix(
        values=[[0.9, 0.2], [0.1, 0.8]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[[0.7, 0.4], [0.3, 0.6]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)

    assert interpreted.execution_config.substrate_support_cutoff == pytest.approx(0.5)
    assert interpreted.execution_config.network_correlation_threshold == pytest.approx(
        0.5
    )
    assert interpreted.execution_config.network_policy == "signed"
    assert (
        interpreted.execution_config.assignment_policy
        == SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
    )
    assert interpreted.execution_config.score_preconditioning_policy == (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )
    assert interpreted.execution_config.requested_module_count is None
    assert interpreted.execution_config.module_selection_primary_threshold == (
        pytest.approx(0.5)
    )
    assert interpreted.execution_config.module_selection_fallback_threshold == (
        pytest.approx(0.1)
    )
    assert interpreted.execution_config.module_selection_max_clusters == 10
    assert interpreted.execution_config.cluster_tree_backend == (
        SIGNALOME_CLUSTER_TREE_BACKEND_EXACT
    )
    assert interpreted.execution_config.candidate_scoring_backend == (
        SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL
    )
    assert interpreted.execution_config.max_exact_cluster_tree_sites == (
        SIGNALOME_MAX_EXACT_CLUSTER_TREE_SITES_DEFAULT
    )
    assert interpreted.execution_config.max_full_correlation_sites == (
        SIGNALOME_MAX_FULL_CORRELATION_SITES_DEFAULT
    )
    assert interpreted.execution_config.clustering_backend_alias == (
        SIGNALOME_CLUSTERING_BACKEND_EXACT
    )


def test_interpreter_preconditions_downstream_scores_without_dropping_prediction_rows() -> (
    None
):
    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1],
            [0.2, 0.8],
            [0.7, 0.3],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [float("nan"), float("nan")],
            [1.0, float("nan")],
            [2.0, 3.0],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)
    assert interpreted.prediction_matrix.index.tolist() == site_ids
    assert interpreted.downstream_score_matrix.index.tolist() == [
        "P2;S2;",
        "P3;S3;",
    ]
    assert pd.isna(interpreted.downstream_score_matrix.loc["P2;S2;", "K2"])
    assert interpreted.score_preconditioning_diagnostics.input_row_count == 3
    assert (
        interpreted.score_preconditioning_diagnostics.dropped_all_missing_row_count == 1
    )
    assert interpreted.score_preconditioning_diagnostics.retained_row_count == 2
    assert interpreted.score_preconditioning_diagnostics.policy == (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )


def test_interpreter_reports_zero_drop_preconditioning_diagnostics() -> None:
    site_ids = ["P1;S1;", "P2;S2;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1],
            [0.2, 0.8],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, float("nan")],
            [2.0, 3.0],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)
    assert interpreted.downstream_score_matrix.index.tolist() == site_ids
    assert interpreted.score_preconditioning_diagnostics.input_row_count == 2
    assert (
        interpreted.score_preconditioning_diagnostics.dropped_all_missing_row_count == 0
    )
    assert interpreted.score_preconditioning_diagnostics.retained_row_count == 2
    assert interpreted.score_preconditioning_diagnostics.policy == (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )


def test_interpreter_respects_explicit_allow_and_report_preconditioning_policy() -> (
    None
):
    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1],
            [0.2, 0.8],
            [0.7, 0.3],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [float("nan"), float("nan")],
            [1.0, float("nan")],
            [2.0, 3.0],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(
            substrate_support_cutoff=0.5,
            score_preconditioning_policy=(
                SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
            ),
        ),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)
    assert interpreted.downstream_score_matrix.index.tolist() == [
        "P2;S2;",
        "P3;S3;",
    ]
    assert interpreted.score_preconditioning_diagnostics.input_row_count == 3
    assert (
        interpreted.score_preconditioning_diagnostics.dropped_all_missing_row_count == 1
    )
    assert interpreted.score_preconditioning_diagnostics.retained_row_count == 2
    assert interpreted.score_preconditioning_diagnostics.policy == (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )


def test_interpreter_fails_when_error_on_drop_policy_detects_all_missing_rows() -> None:
    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1],
            [0.2, 0.8],
            [0.7, 0.3],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [float("nan"), float("nan")],
            [1.0, float("nan")],
            [2.0, 3.0],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(
            substrate_support_cutoff=0.5,
            score_preconditioning_policy=(
                SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
            ),
        ),
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowInterpreter().run(request)

    error = exc_info.value
    message = str(error)
    assert error.seam == SIGNALOME_INTERPRETER_SCORE_PRECONDITIONING_SEAM
    assert error.details["aligned_score_sites"] == 3
    assert error.details["aligned_score_kinases"] == 2
    assert error.details["dropped_all_missing_row_count"] == 1
    assert error.details["retained_row_count"] == 2
    assert error.details["score_preconditioning_policy"] == (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
    )
    assert "dropped_all_missing_row_count=1" in message
    assert (
        "score_preconditioning_policy="
        f"{SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP}"
    ) in message


def test_interpreter_allows_error_on_drop_policy_when_no_rows_require_drop() -> None:
    site_ids = ["P1;S1;", "P2;S2;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1],
            [0.2, 0.8],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, float("nan")],
            [2.0, 3.0],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(
            substrate_support_cutoff=0.5,
            score_preconditioning_policy=(
                SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
            ),
        ),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)
    assert interpreted.downstream_score_matrix.index.tolist() == site_ids
    assert interpreted.score_preconditioning_diagnostics.input_row_count == 2
    assert (
        interpreted.score_preconditioning_diagnostics.dropped_all_missing_row_count == 0
    )
    assert interpreted.score_preconditioning_diagnostics.retained_row_count == 2
    assert interpreted.score_preconditioning_diagnostics.policy == (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
    )


def test_executor_uses_preconditioned_scores_when_missing_rows_are_present() -> None:
    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1],
            [0.2, 0.8],
            [0.7, 0.3],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [float("nan"), float("nan")],
            [1.0, 2.0],
            [2.0, 4.0],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)
    result = SignalomeWorkflowExecutor().run(interpreted)
    assert not result.kinase_network.edges.empty
    assert result.score_preconditioning_diagnostics.input_row_count == 3
    assert result.score_preconditioning_diagnostics.dropped_all_missing_row_count == 1
    assert result.score_preconditioning_diagnostics.retained_row_count == 2
    assert result.score_preconditioning_diagnostics.policy == (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )
    assert result.provenance is not None
    scale_guard = result.provenance.workflow_parameters["scale_guard"]
    assert scale_guard == {
        "site_count": 2,
        "cluster_tree_backend": SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
        "candidate_scoring_backend": SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
        "max_exact_cluster_tree_sites": SIGNALOME_MAX_EXACT_CLUSTER_TREE_SITES_DEFAULT,
        "max_full_correlation_sites": SIGNALOME_MAX_FULL_CORRELATION_SITES_DEFAULT,
        "exact_cluster_tree_built": True,
        "candidate_scoring_mode": SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
        "candidate_scoring_sampling": None,
        "scale_guard_passed": True,
    }


def test_sampled_candidate_scoring_records_sampling_provenance() -> None:
    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1],
            [0.1, 0.9],
            [0.8, 0.8],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, 0.0],
            [0.0, 1.0],
            [0.8, 0.7],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(
            substrate_support_cutoff=0.5,
            cluster_tree_backend=SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
            candidate_scoring_backend=SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
        ),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)
    result = SignalomeWorkflowExecutor().run(interpreted)

    assert result.provenance is not None
    scale_guard = result.provenance.workflow_parameters["scale_guard"]
    assert (
        scale_guard["candidate_scoring_backend"]
        == SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED
    )
    sampled = scale_guard["candidate_scoring_sampling"]
    assert isinstance(sampled, dict)
    assert sampled["sampling_cap"] == MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER
    assert sampled["sampling_method"] == SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD
    assert (
        sampled["deterministic_seed_policy"]
        == SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY
    )
    assert int(sampled["actual_sampled_pair_count"]) >= 0
    per_cluster_summary = sampled["per_cluster_sample_count_summary"]
    assert isinstance(per_cluster_summary, dict)
    assert {"min", "max", "mean", "total"} <= set(per_cluster_summary)
    assert int(per_cluster_summary["min"]) >= 0
    assert int(per_cluster_summary["max"]) >= int(per_cluster_summary["min"])
    assert float(per_cluster_summary["mean"]) >= 0.0
    assert int(per_cluster_summary["total"]) >= int(per_cluster_summary["max"])


def test_sampled_backend_records_zero_sampling_provenance_when_scoring_not_evaluated() -> (
    None
):
    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1],
            [0.1, 0.9],
            [0.8, 0.8],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, 0.0],
            [0.0, 1.0],
            [0.8, 0.7],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(
            substrate_support_cutoff=0.5,
            module_count=2,
            cluster_tree_backend=SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
            candidate_scoring_backend=SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
        ),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)
    result = SignalomeWorkflowExecutor().run(interpreted)

    assert result.provenance is not None
    sampled = result.provenance.workflow_parameters["scale_guard"][
        "candidate_scoring_sampling"
    ]
    assert sampled == {
        "sampling_cap": MAX_APPROX_CORRELATION_SAMPLES_PER_CLUSTER,
        "sampling_method": SIGNALOME_CANDIDATE_SCORING_SAMPLING_METHOD,
        "deterministic_seed_policy": (SIGNALOME_CANDIDATE_SCORING_SAMPLING_SEED_POLICY),
        "actual_sampled_pair_count": 0,
        "per_cluster_sample_count_summary": {
            "min": 0,
            "max": 0,
            "mean": 0.0,
            "total": 0,
        },
    }


def test_exact_cluster_tree_guard_fails_before_build_cluster_tree_call(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.signalomes.clustering as clustering_module

    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.95, 0.1],
            [0.1, 0.95],
            [0.8, 0.7],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, 0.0],
            [0.0, 1.0],
            [0.8, 0.7],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(
            substrate_support_cutoff=0.5,
            module_count=2,
            cluster_tree_backend=SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
            candidate_scoring_backend=SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
            max_exact_cluster_tree_sites=2,
        ),
    )
    interpreted = SignalomeWorkflowInterpreter().run(request)

    tree_calls: list[str] = []

    def _build_tree_should_not_run(scoring_values: object) -> object:
        del scoring_values
        tree_calls.append("called")
        raise AssertionError(
            "build_cluster_tree should not run when exact tree guard fails"
        )

    monkeypatch.setattr(
        clustering_module,
        "build_cluster_tree",
        _build_tree_should_not_run,
    )

    with pytest.raises(SignalomeScaleError) as exc_info:
        SignalomeWorkflowExecutor().run(interpreted)

    message = str(exc_info.value).lower()
    assert "3 sites" in message
    assert "max_exact_cluster_tree_sites=2" in message
    assert "cluster_tree_backend='exact'" in message
    assert "exact cluster-tree construction" in message
    assert tree_calls == []


def test_sampled_candidate_scoring_still_hits_exact_tree_guard_without_approx_tree(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.signalomes.clustering as clustering_module

    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.95, 0.1],
            [0.1, 0.95],
            [0.8, 0.7],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, 0.0],
            [0.0, 1.0],
            [0.8, 0.7],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(
            substrate_support_cutoff=0.5,
            module_count=2,
            cluster_tree_backend=SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
            candidate_scoring_backend=SIGNALOME_CANDIDATE_SCORING_BACKEND_SAMPLED,
            max_exact_cluster_tree_sites=1,
        ),
    )
    interpreted = SignalomeWorkflowInterpreter().run(request)

    tree_calls: list[str] = []

    def _build_tree_should_not_run(scoring_values: object) -> object:
        del scoring_values
        tree_calls.append("called")
        raise AssertionError(
            "build_cluster_tree should not run when exact tree guard fails"
        )

    monkeypatch.setattr(
        clustering_module,
        "build_cluster_tree",
        _build_tree_should_not_run,
    )

    with pytest.raises(SignalomeScaleError) as exc_info:
        SignalomeWorkflowExecutor().run(interpreted)

    message = str(exc_info.value).lower()
    assert "max_exact_cluster_tree_sites=1" in message
    assert "candidate_scoring_backend='sampled'" in message
    assert tree_calls == []


def test_deprecated_approximate_alias_is_still_guarded_by_exact_tree_limit(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.signalomes.clustering as clustering_module

    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.95, 0.1],
            [0.1, 0.95],
            [0.8, 0.7],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, 0.0],
            [0.0, 1.0],
            [0.8, 0.7],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    tree_calls: list[str] = []

    def _build_tree_should_not_run(scoring_values: object) -> object:
        del scoring_values
        tree_calls.append("called")
        raise AssertionError(
            "build_cluster_tree should not run when exact tree guard fails"
        )

    monkeypatch.setattr(
        clustering_module,
        "build_cluster_tree",
        _build_tree_should_not_run,
    )

    with warnings.catch_warnings(record=True) as deprecated_warnings:
        warnings.simplefilter("always", DeprecationWarning)
        request = SignalomeWorkflowRequest(
            kinase_result=_kinase_result(
                dataset=dataset,
                prediction_matrix=prediction_matrix,
                score_matrix=score_matrix,
            ),
            config=SignalomeConfig(
                substrate_support_cutoff=0.5,
                module_count=2,
                clustering_backend=SIGNALOME_CLUSTERING_BACKEND_APPROXIMATE,
                max_exact_clustering_sites=1,
            ),
        )
        interpreted = SignalomeWorkflowInterpreter().run(request)
        with pytest.raises(SignalomeScaleError) as exc_info:
            SignalomeWorkflowExecutor().run(interpreted)

    warning_messages = [
        str(warning.message)
        for warning in deprecated_warnings
        if issubclass(warning.category, DeprecationWarning)
    ]
    assert any(
        "config.clustering_backend is deprecated" in msg for msg in warning_messages
    )
    assert any(
        "config.max_exact_clustering_sites is deprecated" in msg
        for msg in warning_messages
    )

    message = str(exc_info.value).lower()
    assert "max_exact_cluster_tree_sites=1" in message
    assert "cluster_tree_backend='exact'" in message
    assert "candidate_scoring_backend='sampled'" in message
    assert tree_calls == []


def test_full_candidate_scoring_guard_triggers_above_full_correlation_limit() -> None:
    site_ids = ["P1;S1;", "P2;S2;", "P3;S3;"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[
            [0.95, 0.1],
            [0.1, 0.95],
            [0.8, 0.7],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, 0.0],
            [0.0, 1.0],
            [0.8, 0.7],
        ],
        site_ids=site_ids,
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(
            substrate_support_cutoff=0.5,
            cluster_tree_backend=SIGNALOME_CLUSTER_TREE_BACKEND_EXACT,
            candidate_scoring_backend=SIGNALOME_CANDIDATE_SCORING_BACKEND_FULL,
            max_exact_cluster_tree_sites=10,
            max_full_correlation_sites=2,
        ),
    )
    interpreted = SignalomeWorkflowInterpreter().run(request)

    with pytest.raises(SignalomeScaleError) as exc_info:
        SignalomeWorkflowExecutor().run(interpreted)

    message = str(exc_info.value).lower()
    assert "max_full_correlation_sites=2" in message
    assert "candidate_scoring_backend='sampled'" in message


def test_signalome_grouping_does_not_collapse_distinct_protein_ids_with_shared_gene_symbol() -> (
    None
):
    dataset = _dataset(
        site_ids=["MAPK14;S1;", "MAPK14;S2;"],
        gene_symbols=["MAPK14", "MAPK14"],
        protein_ids=["P28482-1", "P28482-2"],
    )
    prediction_matrix = _matrix(
        values=[[0.9, 0.1], [0.1, 0.9]],
        site_ids=["MAPK14;S1;", "MAPK14;S2;"],
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[[1.0, 0.0], [0.0, 1.0]],
        site_ids=["MAPK14;S1;", "MAPK14;S2;"],
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    resolved = SignalomeWorkflowInterpreter().run(request)
    result = SignalomeWorkflowExecutor().run(resolved)
    assignments = result.module_assignments.table
    proteins = assignments.loc[:, "protein_id"].tolist()
    assert proteins == ["P28482-1", "P28482-2"]
    assert assignments.loc[:, "module_id"].astype("int64").ge(0).all()


def test_boundary_error_reports_no_support_cutoff_support_counts() -> None:
    dataset = _dataset(site_ids=["P1;S1;", "P2;S2;"])
    prediction_matrix = _matrix(
        values=[[0.2, 0.4], [0.3, 0.1]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[[1.0, 2.0], [3.0, 4.0]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.9),
    )
    interpreted = SignalomeWorkflowInterpreter().run(request)

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowExecutor().run(interpreted)

    message = str(exc_info.value)
    assert f"seam={SIGNALOME_EXECUTOR_KINASE_SUPPORT_SEAM}" in message
    assert "prediction_sites=2" in message
    assert "prediction_kinases=2" in message
    assert "supported_sites=0" in message
    assert "supported_kinases=0" in message
    assert "substrate_support_cutoff=0.9" in message


def test_boundary_error_reports_module_construction_degeneracy_counts() -> None:
    dataset = _dataset(site_ids=["P1;S1;", "P2;S2;"])
    prediction_matrix = _matrix(
        values=[[0.9], [0.8]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1"],
    )
    score_matrix = _matrix(
        values=[[1.0], [2.0]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )
    interpreted = SignalomeWorkflowInterpreter().run(request)

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowExecutor().run(interpreted)

    message = str(exc_info.value)
    assert f"seam={SIGNALOME_EXECUTOR_MODULE_CONSTRUCTION_SEAM}" in message
    assert "module_count=1" in message
    assert "supported_kinases=1" in message
    assert "prediction_kinases=1" in message
    assert "substrate_support_cutoff=0.5" in message
    assert "network_correlation_threshold=0.5" in message


def test_boundary_error_reports_network_failure_modes() -> None:
    dataset = _dataset(site_ids=["P1;S1;", "P2;S2;"])
    prediction_matrix = _matrix(
        values=[[0.9, 0.1], [0.1, 0.9]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1", "K2"],
    )

    score_matrix_missing_kinase = _matrix(
        values=[[1.0], [2.0]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1"],
    )
    resolved_missing_kinase = ResolvedSignalomeWorkflowRequest(
        dataset=dataset,
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix_missing_kinase,
        ),
        execution_config=_execution_config(
            SignalomeConfig(substrate_support_cutoff=0.5)
        ),
        downstream_score_matrix=score_matrix_missing_kinase,
        downstream_score_source="rank_weighted_fusion_scores",
        prediction_matrix=prediction_matrix,
        site_to_protein=pd.Series(
            ["P1", "P2"],
            index=pd.Index(["P1;S1;", "P2;S2;"], name="site_id"),
            name="protein_id",
            dtype=str,
        ),
    )

    with pytest.raises(WorkflowBoundaryError) as missing_exc:
        SignalomeWorkflowExecutor().run(resolved_missing_kinase)

    missing_message = str(missing_exc.value)
    assert f"seam={SIGNALOME_EXECUTOR_NETWORK_SEAM}" in missing_message
    assert "shared_kinases=2" in missing_message
    assert "supported_kinases=2" in missing_message
    assert (
        "stage_error=downstream score matrix is missing kinases required for signalome network"
        in missing_message
    )

    score_matrix_zero_variance = _matrix(
        values=[[1.0, 2.0], [1.0, 2.0]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix_zero_variance,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )
    interpreted = SignalomeWorkflowInterpreter().run(request)

    with pytest.raises(WorkflowBoundaryError) as variance_exc:
        SignalomeWorkflowExecutor().run(interpreted)

    variance_message = str(variance_exc.value)
    assert f"seam={SIGNALOME_EXECUTOR_NETWORK_SEAM}" in variance_message
    assert "shared_kinases=2" in variance_message
    assert "supported_kinases=2" in variance_message
    assert "downstream_score_sites=2" in variance_message
    assert "score_variance_kinases=0" in variance_message
    assert "network_correlation_threshold=0.5" in variance_message


def test_boundary_error_reports_expanded_signalome_failure_seam(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.workflows.signalome.executor as executor_module

    dataset = _dataset(site_ids=["P1;S1;", "P2;S2;", "P3;S3;"])
    prediction_matrix = _matrix(
        values=[
            [0.95, 0.1],
            [0.1, 0.95],
            [0.8, 0.7],
        ],
        site_ids=["P1;S1;", "P2;S2;", "P3;S3;"],
        kinases=["K1", "K2"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, 0.0],
            [0.0, 1.0],
            [0.8, 0.7],
        ],
        site_ids=["P1;S1;", "P2;S2;", "P3;S3;"],
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )
    interpreted = SignalomeWorkflowInterpreter().run(request)

    def _raise_expanded_stage_error(**_: object) -> pd.DataFrame:
        raise WorkflowStageError("expanded signalome seam regression test")

    monkeypatch.setattr(
        executor_module,
        "build_expanded_signalome_table",
        _raise_expanded_stage_error,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowExecutor().run(interpreted)

    message = str(exc_info.value)
    assert f"seam={SIGNALOME_EXECUTOR_EXPANDED_SIGNALOME_SEAM}" in message
    assert f"assignment_policy={SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY}" in message
    assert "stage_error=expanded signalome seam regression test" in message


def test_boundary_error_reports_context_table_site_membership_failure_seam(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.workflows.signalome.executor as executor_module

    interpreted = _interpreted_request_for_context_failure()

    def _raise_site_membership_error(**_: object) -> pd.DataFrame:
        raise ValueError("site membership context seam regression test")

    monkeypatch.setattr(
        executor_module,
        "build_site_membership_table",
        _raise_site_membership_error,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowExecutor().run(interpreted)

    error = exc_info.value
    message = str(error)
    assert error.seam == SIGNALOME_EXECUTOR_CONTEXT_TABLES_SEAM
    assert "prediction_sites=3" in message
    assert "prediction_kinases=2" in message
    assert "module_assignment_rows=3" in message
    assert "module_assignment_columns=" in message
    assert "site_cluster_count=" in message
    assert "supported_sites=3" in message
    assert "supported_kinases=2" in message
    assert "substrate_support_cutoff=0.5" in message
    assert f"assignment_policy={SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY}" in message
    assert "stage_error=site membership context seam regression test" in message


def test_boundary_error_reports_context_table_protein_context_failure_seam(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.workflows.signalome.executor as executor_module

    interpreted = _interpreted_request_for_context_failure()

    def _raise_protein_context_error(**_: object) -> pd.DataFrame:
        raise WorkflowStageError("protein context seam regression test")

    monkeypatch.setattr(
        executor_module,
        "build_protein_site_context_table",
        _raise_protein_context_error,
    )

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowExecutor().run(interpreted)

    error = exc_info.value
    message = str(error)
    assert error.seam == SIGNALOME_EXECUTOR_CONTEXT_TABLES_SEAM
    assert "prediction_sites=3" in message
    assert "prediction_kinases=2" in message
    assert "module_assignment_rows=3" in message
    assert "module_assignment_columns=" in message
    assert "site_cluster_count=" in message
    assert "supported_sites=3" in message
    assert "supported_kinases=2" in message
    assert "substrate_support_cutoff=0.5" in message
    assert f"assignment_policy={SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY}" in message
    assert "stage_error=protein context seam regression test" in message


def test_signalome_result_rejects_malformed_site_membership_immediately() -> None:
    dataset = _dataset(site_ids=["P1;S1;"])
    prediction_matrix = _matrix(
        values=[[0.9]],
        site_ids=["P1;S1;"],
        kinases=["K1"],
    )
    score_matrix = _matrix(
        values=[[0.8]],
        site_ids=["P1;S1;"],
        kinases=["K1"],
    )
    kinase_result = _kinase_result(
        dataset=dataset,
        prediction_matrix=prediction_matrix,
        score_matrix=score_matrix,
    )
    with pytest.raises(WorkflowValidationError, match="missing required columns"):
        SignalomeWorkflowResult(
            dataset=dataset,
            kinase_result=kinase_result,
            module_assignments=SignalomeAssignments(
                table=pd.DataFrame(
                    {"protein_id": ["P1"], "module_id": [1]},
                    index=pd.Index(["P1;S1;"], name="site_id"),
                )
            ),
            signalome_modules=SignalomeModules(
                table=pd.DataFrame(
                    {"K1": [100.0]},
                    index=pd.Index([1], name="module_id"),
                )
            ),
            kinase_network=KinaseNetwork(
                edges=pd.DataFrame(
                    columns=["source_kinase", "target_kinase", "correlation"]
                )
            ),
            site_membership=pd.DataFrame({"site_id": ["P1;S1;"]}),
        )


def test_support_cutoff_changes_substrate_support_without_changing_network_edges() -> (
    None
):
    dataset = _dataset(site_ids=["P1;S1;", "P1;S2;", "P2;S3;", "P2;S4;"])
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1, 0.6],
            [0.8, 0.2, 0.55],
            [0.2, 0.85, 0.4],
            [0.1, 0.9, 0.35],
        ],
        site_ids=["P1;S1;", "P1;S2;", "P2;S3;", "P2;S4;"],
        kinases=["K1", "K2", "K3"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, 1.0, 4.0],
            [2.0, 2.1, 3.0],
            [3.0, 2.9, 2.0],
            [4.0, 4.1, 1.0],
        ],
        site_ids=["P1;S1;", "P1;S2;", "P2;S3;", "P2;S4;"],
        kinases=["K1", "K2", "K3"],
    )
    kinase_result = _kinase_result(
        dataset=dataset,
        prediction_matrix=prediction_matrix,
        score_matrix=score_matrix,
    )
    executor = SignalomeWorkflowExecutor()

    low_support_resolved = SignalomeWorkflowInterpreter().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(
                substrate_support_cutoff=0.5,
                network_correlation_threshold=0.95,
            ),
        )
    )
    high_support_resolved = SignalomeWorkflowInterpreter().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(
                substrate_support_cutoff=0.75,
                network_correlation_threshold=0.95,
            ),
        )
    )

    low_support = executor.run(low_support_resolved)
    high_support = executor.run(high_support_resolved)

    pd.testing.assert_frame_equal(
        low_support.kinase_network.edges,
        high_support.kinase_network.edges,
        check_dtype=False,
    )
    assert (
        low_support.kinase_network.nodes.loc["K3", "n_substrates"]
        > high_support.kinase_network.nodes.loc["K3", "n_substrates"]
    )


def test_network_threshold_changes_edge_sparsity_without_changing_substrate_support() -> (
    None
):
    dataset = _dataset(site_ids=["P1;S1;", "P1;S2;", "P2;S3;", "P2;S4;"])
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1, 0.6],
            [0.8, 0.2, 0.55],
            [0.2, 0.85, 0.4],
            [0.1, 0.9, 0.35],
        ],
        site_ids=["P1;S1;", "P1;S2;", "P2;S3;", "P2;S4;"],
        kinases=["K1", "K2", "K3"],
    )
    score_matrix = _matrix(
        values=[
            [1.0, 1.0, 4.0],
            [2.0, 2.1, 3.0],
            [3.0, 2.9, 2.0],
            [4.0, 4.1, 1.0],
        ],
        site_ids=["P1;S1;", "P1;S2;", "P2;S3;", "P2;S4;"],
        kinases=["K1", "K2", "K3"],
    )
    kinase_result = _kinase_result(
        dataset=dataset,
        prediction_matrix=prediction_matrix,
        score_matrix=score_matrix,
    )
    executor = SignalomeWorkflowExecutor()

    low_threshold_resolved = SignalomeWorkflowInterpreter().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(
                substrate_support_cutoff=0.5,
                network_correlation_threshold=0.95,
            ),
        )
    )
    high_threshold_resolved = SignalomeWorkflowInterpreter().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(
                substrate_support_cutoff=0.5,
                network_correlation_threshold=0.999,
            ),
        )
    )

    low_threshold = executor.run(low_threshold_resolved)
    high_threshold = executor.run(high_threshold_resolved)

    pd.testing.assert_frame_equal(
        low_threshold.signalome_modules.table,
        high_threshold.signalome_modules.table,
        check_dtype=False,
    )
    pd.testing.assert_frame_equal(
        low_threshold.kinase_network.nodes.loc[:, ["n_substrates"]],
        high_threshold.kinase_network.nodes.loc[:, ["n_substrates"]],
        check_dtype=False,
    )
    assert (
        low_threshold.kinase_network.edges.shape[0]
        > high_threshold.kinase_network.edges.shape[0]
    )


def test_executor_orchestrates_signalome_domain_services(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import phospy.workflows.signalome.executor as executor_module
    from phospy.signalomes.clustering import ClusterSitesResult
    from phospy.signalomes.models import SignalomeModuleSelectionDiagnostics

    site_ids = ["P1;S1;", "P2;S2;"]
    kinases = ["K1", "K2"]
    dataset = _dataset(site_ids=site_ids)
    prediction_matrix = _matrix(
        values=[[0.9, 0.1], [0.1, 0.9]],
        site_ids=site_ids,
        kinases=kinases,
    )
    score_matrix = _matrix(
        values=[[1.0, 0.5], [0.2, 1.1]],
        site_ids=site_ids,
        kinases=kinases,
    )
    resolved = ResolvedSignalomeWorkflowRequest(
        dataset=dataset,
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        execution_config=_execution_config(
            SignalomeConfig(substrate_support_cutoff=0.5)
        ),
        downstream_score_matrix=score_matrix,
        downstream_score_source="rank_weighted_fusion_scores",
        prediction_matrix=prediction_matrix,
        site_to_protein=pd.Series(
            ["P1", "P2"],
            index=pd.Index(site_ids, name="site_id"),
            name="protein_id",
            dtype=str,
        ),
    )

    call_order: list[str] = []

    def _cluster(**_: object) -> ClusterSitesResult:
        call_order.append("cluster")
        return ClusterSitesResult(
            site_clusters=pd.Series(
                [1, 2],
                index=prediction_matrix.index.copy(),
                dtype=int,
                name=SITE_CLUSTER_COLUMN,
            ),
            module_selection_diagnostics=SignalomeModuleSelectionDiagnostics(
                strategy="correlation_thresholds",
                selected_module_count=2,
                requested_module_count=None,
                threshold_used=0.5,
                max_clusters_evaluated=2,
                candidate_scores={},
                reason="test stub",
            ),
        )

    def _derive(**_: object) -> pd.Series:
        call_order.append("derive")
        return pd.Series({"P1": 1, "P2": 2}, dtype="int64", name="module_id")

    def _assignments(**_: object) -> pd.DataFrame:
        call_order.append("assignments")
        return pd.DataFrame(
            {
                "protein_id": ["P1", "P2"],
                "module_id": [1, 2],
            },
            index=pd.Index(site_ids, name="site_id"),
        )

    def _substrates(**_: object) -> dict[str, tuple[str, ...]]:
        call_order.append("substrates")
        return {"K1": ("P1;S1;",), "K2": ("P2;S2;",)}

    def _modules(**_: object) -> pd.DataFrame:
        call_order.append("modules")
        return pd.DataFrame(
            {"K1": [100.0, 0.0], "K2": [0.0, 100.0]},
            index=pd.Index([1, 2], name="module_id"),
            dtype=float,
        )

    def _network(
        **_: object,
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, object]:
        call_order.append("network")
        return (
            pd.DataFrame(
                {
                    "source_kinase": ["K1"],
                    "target_kinase": ["K2"],
                    "correlation": [0.95],
                }
            ),
            pd.DataFrame(
                {"degree": [1, 1], "n_substrates": [1, 1]},
                index=pd.Index(kinases, name="kinase"),
            ).astype({"degree": "int64", "n_substrates": "int64"}),
            pd.DataFrame(
                {
                    "source_kinase": ["K1"],
                    "target_kinase": ["K2"],
                    "correlation": [0.95],
                    "correlation_status": ["finite"],
                    "valid_observations": [2],
                    "correlation_reason": [None],
                }
            ),
            executor_module.SignalomeNetworkCorrelationDiagnostics(
                total_candidate_correlations=1,
                finite_correlations=1,
                undefined_correlations=0,
                constant_profile_correlations=0,
                insufficient_observation_correlations=0,
                missing_value_correlations=0,
                non_finite_value_correlations=0,
                edges_created=1,
                edges_skipped_non_finite_correlation=0,
            ),
        )

    def _expanded(**_: object) -> pd.DataFrame:
        call_order.append("expanded")
        return pd.DataFrame(
            {
                "kinase": ["K1"],
                "row_kind": [EXPANDED_SIGNALOME_ROW_KIND_SUMMARY],
                "assignment_policy": [SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY],
                "linked_kinases": ['["K1","K2"]'],
                "regulated_module_ids": ["[1]"],
                "site_id": [""],
                "site_order": [-1],
                "protein_id": [""],
                "module_id": [0],
                "support_kinases": ["[]"],
                "support_weight": [0.0],
                "top_kinase": [""],
                "top_score": [float("nan")],
            }
        )

    monkeypatch.setattr(executor_module, "cluster_sites_with_diagnostics", _cluster)
    monkeypatch.setattr(executor_module, "derive_protein_modules", _derive)
    monkeypatch.setattr(executor_module, "build_module_assignments", _assignments)
    monkeypatch.setattr(executor_module, "select_kinase_substrates", _substrates)
    monkeypatch.setattr(executor_module, "build_signalome_module_table", _modules)
    monkeypatch.setattr(
        executor_module,
        "build_kinase_network_with_diagnostics",
        _network,
    )
    monkeypatch.setattr(executor_module, "build_expanded_signalome_table", _expanded)

    result = SignalomeWorkflowExecutor().run(resolved)

    assert call_order == [
        "cluster",
        "derive",
        "assignments",
        "substrates",
        "modules",
        "network",
        "expanded",
    ]
    assert not result.kinase_network.edges.empty

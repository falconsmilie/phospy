from __future__ import annotations

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyPhosphoDataset,
    KinasePredictionResult,
    KinaseScoringResult,
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    SignalomeConfig,
    SignalomeWorkflowRequest,
)
from phospy.errors import WorkflowBoundaryError
from phospy.transformations.models import TransformationState
from phospy.workflows.signalome.contracts import ResolvedSignalomeWorkflowRequest
from phospy.workflows.signalome.executor import SignalomeWorkflowExecutor
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter


def _dataset(
    *,
    site_ids: list[str],
    gene_symbols: list[str] | None = None,
    protein_ids: list[str] | None = None,
) -> AnalysisReadyPhosphoDataset:
    if gene_symbols is None:
        gene_symbols = [str(site_id).split(";", 1)[0] for site_id in site_ids]
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
        },
        index=site_ids,
    )
    if protein_ids is not None:
        site_metadata.loc[:, "protein_id"] = protein_ids
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        transformation_state=TransformationState.raw(has_total_matrix=False),
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
            combined_scores=(
                score_matrix if combined_score_matrix is None else combined_score_matrix
            ),
        ),
        prediction_result=KinasePredictionResult(pred_mat=prediction_matrix),
        activity_result=None,
    )


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

    message = str(exc_info.value)
    assert "seam=signalome.interpreter.site_alignment" in message
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
    assert "seam=signalome.interpreter.kinase_overlap" in message
    assert "prediction_kinases=2" in message
    assert "score_kinases=2" in message
    assert "shared_kinases=0" in message
    assert "next_action=" in message


def test_boundary_error_reports_unusable_protein_mapping_counts() -> None:
    dataset = _dataset(
        site_ids=[";S1;", ";S2;"],
        gene_symbols=["MAPK14", "GSK3B"],
    )
    prediction_matrix = _matrix(
        values=[[0.9], [0.8]],
        site_ids=[";S1;", ";S2;"],
        kinases=["K1"],
    )
    score_matrix = _matrix(
        values=[[1.0], [2.0]],
        site_ids=[";S1;", ";S2;"],
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
    assert "seam=signalome.interpreter.protein_mapping" in message
    assert "protein_resolution_source=dataset.phospho.index_protein_prefix" in message
    assert "interpreted_sites=2" in message
    assert "resolved_protein_sites=0" in message
    assert "unresolved_protein_sites=2" in message
    assert "next_action=" in message


def test_interpreter_uses_explicit_site_metadata_protein_id_when_present() -> None:
    dataset = _dataset(
        site_ids=[";S1;", ";S2;"],
        gene_symbols=["MAPK14", "MAPK14"],
        protein_ids=["P28482-1", "P28482-2"],
    )
    prediction_matrix = _matrix(
        values=[[0.9], [0.8]],
        site_ids=[";S1;", ";S2;"],
        kinases=["K1"],
    )
    score_matrix = _matrix(
        values=[[1.0], [2.0]],
        site_ids=[";S1;", ";S2;"],
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


def test_interpreter_prefers_combined_scores_for_downstream_signalome_matrix() -> None:
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
    combined_scores = _matrix(
        values=[[0.7, 0.4], [0.3, 0.6]],
        site_ids=["P1;S1;", "P2;S2;"],
        kinases=["K1", "K2"],
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=profile_scores,
            combined_score_matrix=combined_scores,
        ),
        config=SignalomeConfig(substrate_support_cutoff=0.5),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)
    pd.testing.assert_frame_equal(
        interpreted.downstream_score_matrix,
        combined_scores,
        check_dtype=False,
    )
    assert interpreted.downstream_score_source == "combined_scores"


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
    assert assignments.loc[:, "module_id"].nunique() == 2


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
    assert "seam=signalome.executor.kinase_support" in message
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
    assert "seam=signalome.executor.module_construction" in message
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
        config=SignalomeConfig(substrate_support_cutoff=0.5),
        downstream_score_matrix=score_matrix_missing_kinase,
        downstream_score_source="combined_scores",
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
    assert "seam=signalome.executor.network" in missing_message
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
    assert "seam=signalome.executor.network" in variance_message
    assert "shared_kinases=2" in variance_message
    assert "supported_kinases=2" in variance_message
    assert "downstream_score_sites=2" in variance_message
    assert "score_variance_kinases=0" in variance_message
    assert "network_correlation_threshold=0.5" in variance_message


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

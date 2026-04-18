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
) -> KinaseWorkflowResult:
    return KinaseWorkflowResult(
        dataset=dataset,
        references=_bundle(site_ids=dataset.phospho.index.astype(str).tolist()),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            combined_scores=score_matrix,
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
        config=SignalomeConfig(signalome_cutoff=0.5),
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
        config=SignalomeConfig(signalome_cutoff=0.5),
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
    dataset = _dataset(site_ids=[";S1;", ";S2;"], gene_symbols=["", " "])
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
        config=SignalomeConfig(signalome_cutoff=0.5),
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
        config=SignalomeConfig(signalome_cutoff=0.5),
    )

    interpreted = SignalomeWorkflowInterpreter().run(request)
    assert interpreted.site_to_protein.tolist() == ["P28482-1", "P28482-2"]


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
        config=SignalomeConfig(signalome_cutoff=0.5),
    )

    resolved = SignalomeWorkflowInterpreter().run(request)
    result = SignalomeWorkflowExecutor().run(resolved)
    assignments = result.module_assignments.table
    proteins = assignments.loc[:, "protein_id"].tolist()
    assert proteins == ["P28482-1", "P28482-2"]
    assert assignments.loc[:, "module_id"].nunique() == 2


def test_boundary_error_reports_no_cutoff_support_counts() -> None:
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
        config=SignalomeConfig(signalome_cutoff=0.9),
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
    assert "signalome_cutoff=0.9" in message


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
        config=SignalomeConfig(signalome_cutoff=0.5),
    )
    interpreted = SignalomeWorkflowInterpreter().run(request)

    with pytest.raises(WorkflowBoundaryError) as exc_info:
        SignalomeWorkflowExecutor().run(interpreted)

    message = str(exc_info.value)
    assert "seam=signalome.executor.module_construction" in message
    assert "module_count=1" in message
    assert "supported_kinases=1" in message
    assert "prediction_kinases=1" in message
    assert "signalome_cutoff=0.5" in message


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
        config=SignalomeConfig(signalome_cutoff=0.5),
        score_matrix=score_matrix_missing_kinase,
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
        "stage_error=score matrix is missing kinases required for signalome network"
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
        config=SignalomeConfig(signalome_cutoff=0.5),
    )
    interpreted = SignalomeWorkflowInterpreter().run(request)

    with pytest.raises(WorkflowBoundaryError) as variance_exc:
        SignalomeWorkflowExecutor().run(interpreted)

    variance_message = str(variance_exc.value)
    assert "seam=signalome.executor.network" in variance_message
    assert "shared_kinases=2" in variance_message
    assert "supported_kinases=2" in variance_message
    assert "score_sites=2" in variance_message
    assert "score_variance_kinases=0" in variance_message
    assert "signalome_cutoff=0.5" in variance_message

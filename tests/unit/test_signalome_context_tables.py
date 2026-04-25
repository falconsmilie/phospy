from __future__ import annotations

import pandas as pd
import pandas.testing as pdt

from phospy import AnalysisReadyPhosphoDataset
from phospy.api import (
    KinaseWorkflowResult,
    Organism,
    ReferenceBundle,
    SignalomeConfig,
    SignalomeWorkflowRequest,
)
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    SignalomeWorkflowResult,
)
from phospy.signalomes.clustering import (
    cluster_sites_with_diagnostics,
    derive_protein_modules,
)
from phospy.signalomes.science import (
    build_module_assignments,
    build_signalome_module_table,
    select_kinase_substrates,
)
from phospy.workflows.signalome.contracts import ResolvedSignalomeWorkflowRequest
from phospy.workflows.signalome.executor import SignalomeWorkflowExecutor
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter
from tests.support.transformation_states import supported_linear_state


def _dataset() -> AnalysisReadyPhosphoDataset:
    site_ids = ["P1;S1;", "P1;S2;", "P2;S3;", "P3;S4;"]
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0, 4.0],
            "sample_b": [1.2, 2.2, 3.2, 4.2],
        },
        index=site_ids,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["P1", "P1", "P2", "P3"],
            "site": ["S1", "S2", "S3", "S4"],
            "site_sequence": ["A" * 31, "B" * 31, "C" * 31, "D" * 31],
            "protein_id": ["P1", "P1", "P2", "P3"],
        },
        index=site_ids,
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        transformation_state=supported_linear_state(has_total_matrix=False),
    )


def _bundle(site_ids: list[str]) -> ReferenceBundle:
    unique_sites = pd.Index([str(site_id) for site_id in site_ids], name="site_id")
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["K1", "K2"],
                "substrate_site": [unique_sites[0], unique_sites[1]],
            }
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


def _run_signalome_executor() -> tuple[
    SignalomeWorkflowResult, ResolvedSignalomeWorkflowRequest
]:
    dataset = _dataset()
    site_ids = dataset.phospho.index.astype(str).tolist()
    kinases = ["K1", "K2"]
    prediction_matrix = _matrix(
        values=[
            [0.9, 0.1],
            [0.1, 0.9],
            [0.8, 0.2],
            [0.7, 0.3],
        ],
        site_ids=site_ids,
        kinases=kinases,
    )
    score_matrix = _matrix(
        values=[
            [1.0, 0.0],
            [0.0, 1.0],
            [0.9, 0.1],
            [float("nan"), float("nan")],
        ],
        site_ids=site_ids,
        kinases=kinases,
    )
    request = SignalomeWorkflowRequest(
        kinase_result=_kinase_result(
            dataset=dataset,
            prediction_matrix=prediction_matrix,
            score_matrix=score_matrix,
        ),
        config=SignalomeConfig(
            substrate_support_cutoff=0.5,
            network_correlation_threshold=0.3,
            module_count=2,
        ),
    )
    interpreted = SignalomeWorkflowInterpreter().run(request)
    result = SignalomeWorkflowExecutor().run(interpreted)
    return result, interpreted


def test_signalome_result_exposes_site_and_protein_context_tables() -> None:
    result, _ = _run_signalome_executor()

    assert result.site_membership is not None
    assert result.protein_site_context is not None
    assert not result.module_assignments.table.empty
    assert {
        "site_id",
        "protein_id",
        "site_cluster",
        "protein_module_id",
        "included_in_module_table",
        "excluded_reason",
    }.issubset(result.site_membership.columns)
    assert {
        "protein_id",
        "n_sites",
        "site_ids",
        "site_clusters",
        "n_distinct_site_clusters",
        "protein_module_id",
        "multi_site_protein",
        "ambiguous_module_context",
    }.issubset(result.protein_site_context.columns)


def test_site_membership_contains_all_interpreted_site_rows() -> None:
    result, interpreted = _run_signalome_executor()
    assert result.site_membership is not None
    site_membership = result.site_membership

    expected_site_ids = interpreted.prediction_matrix.index.astype(str).tolist()
    assert site_membership.shape[0] == len(expected_site_ids)
    assert set(site_membership.loc[:, "site_id"].astype(str)) == set(expected_site_ids)

    dropped_site = site_membership.set_index("site_id").loc["P3;S4;"]
    assert pd.isna(dropped_site["site_cluster"])
    assert not bool(dropped_site["included_in_module_table"])
    assert (
        str(dropped_site["excluded_reason"]) == "dropped_all_missing_downstream_scores"
    )


def test_protein_site_context_flags_multi_site_and_ambiguity() -> None:
    result, _ = _run_signalome_executor()
    assert result.protein_site_context is not None
    context = result.protein_site_context.set_index("protein_id")

    assert bool(context.loc["P1", "multi_site_protein"])
    assert bool(context.loc["P1", "ambiguous_module_context"])
    assert int(context.loc["P1", "n_sites"]) == 2
    assert int(context.loc["P1", "n_distinct_site_clusters"]) >= 2


def test_context_tables_do_not_change_module_outputs() -> None:
    result, interpreted = _run_signalome_executor()
    config = interpreted.execution_config

    clustering_result = cluster_sites_with_diagnostics(
        scoring_matrix=interpreted.downstream_score_matrix,
        requested_module_count=config.requested_module_count,
        primary_threshold=config.module_selection_primary_threshold,
        fallback_threshold=config.module_selection_fallback_threshold,
        max_clusters=config.module_selection_max_clusters,
    )
    protein_modules = derive_protein_modules(
        site_clusters=clustering_result.site_clusters,
        site_to_protein=interpreted.site_to_protein,
    )
    expected_assignments = build_module_assignments(
        prediction_matrix=interpreted.prediction_matrix,
        site_to_protein=interpreted.site_to_protein,
        protein_modules=protein_modules,
    )
    expected_substrates = select_kinase_substrates(
        prediction_matrix=interpreted.prediction_matrix,
        cutoff=config.substrate_support_cutoff,
    )
    expected_modules = build_signalome_module_table(
        module_assignments=expected_assignments,
        kinase_substrates=expected_substrates,
        kinase_order=interpreted.prediction_matrix.columns.astype(str).tolist(),
        assignment_policy=config.assignment_policy,
    )

    pdt.assert_frame_equal(
        result.module_assignments.table,
        expected_assignments,
        check_dtype=False,
    )
    pdt.assert_frame_equal(
        result.signalome_modules.table,
        expected_modules,
        check_dtype=False,
    )

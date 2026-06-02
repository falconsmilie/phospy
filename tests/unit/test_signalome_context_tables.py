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
from phospy.api.results import (
    KinasePredictionResult,
    KinaseScoringResult,
    SignalomeWorkflowResult,
)
from phospy.errors.workflows import WorkflowStageError
from phospy.science.signalomes.clustering import (
    cluster_sites_with_diagnostics,
    derive_protein_modules,
)
from phospy.science.signalomes.context import (
    build_protein_site_context_table,
    build_site_membership_table,
)
from phospy.science.signalomes.science import (
    build_module_assignments,
    build_signalome_module_table,
    select_kinase_substrates,
)
from phospy.science.sites.site_keys import (
    ProteinScopedPhosphositeKey,
    encode_site_key,
)
from phospy.workflows.signalome.contracts import ResolvedSignalomeWorkflowRequest
from phospy.workflows.signalome.executor import SignalomeWorkflowExecutor
from phospy.workflows.signalome.interpreter import SignalomeWorkflowInterpreter
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import site_key_context_columns


def _site_key(*, protein_id: str, site: str) -> str:
    return encode_site_key(
        ProteinScopedPhosphositeKey(
            organism=Organism.RAT.value,
            protein_namespace="protein_id",
            protein_identifier=protein_id,
            residue=site[0],
            position=int(site[1:]),
        )
    )


def _dataset() -> AnalysisReadyPhosphoDataset:
    protein_ids = ["P1", "P1", "P2", "P3"]
    sites = ["S1", "S2", "S3", "S4"]
    display_ids = [
        f"{protein_id};{site};"
        for protein_id, site in zip(protein_ids, sites, strict=True)
    ]
    site_keys = [
        _site_key(protein_id=protein_id, site=site)
        for protein_id, site in zip(protein_ids, sites, strict=True)
    ]
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 3.0, 4.0],
            "sample_b": [1.2, 2.2, 3.2, 4.2],
        },
        index=pd.Index(site_keys, name="site_key"),
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_keys,
            "display_id": display_ids,
            **site_key_context_columns(site_keys),
            "gene_symbol": protein_ids,
            "site": sites,
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15) for site in sites
            ],
            "protein_id": protein_ids,
        },
        index=pd.Index(site_keys, name="site_key"),
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
            {
                "kinase": ["K1", "K2"],
                "substrate_site": [unique_sites[0], unique_sites[1]],
            }
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
    values: list[list[float]],
    site_ids: list[str],
    kinases: list[str],
) -> pd.DataFrame:
    return pd.DataFrame(
        values,
        index=pd.Index(site_ids, name="site_key"),
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
        references=_bundle(
            site_ids=dataset.site_metadata.loc[:, "display_id"].astype(str).tolist()
        ),
        scoring_result=KinaseScoringResult(
            profile_scores=score_matrix,
            rank_weighted_fusion_scores=score_matrix,
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
        config=build_signalome_config(
            substrate_support_cutoff=0.5,
            network_correlation_threshold=0.3,
            module_count=2,
            score_preconditioning_policy="allow_and_report",
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


def test_site_membership_contains_only_retained_interpreted_sites() -> None:
    result, interpreted = _run_signalome_executor()
    assert result.site_membership is not None
    site_membership = result.site_membership

    expected_site_keys = interpreted.prediction_matrix.index.astype(str).tolist()
    assert site_membership.shape[0] == len(expected_site_keys)
    assert set(site_membership.loc[:, "site_key"].astype(str)) == set(
        expected_site_keys
    )

    assert "P3;S4;" not in set(site_membership.loc[:, "site_id"].astype(str))
    assert interpreted.prediction_matrix.index.equals(
        interpreted.downstream_score_matrix.index
    )
    assert interpreted.site_to_protein.index.equals(
        interpreted.downstream_score_matrix.index
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
        site_metadata=interpreted.dataset._borrow_site_metadata_frame(),
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


def test_site_membership_builder_rejects_non_empty_missing_required_columns() -> None:
    with pytest.raises(WorkflowStageError, match="missing columns"):
        build_site_membership_table(
            module_assignments=pd.DataFrame(
                {"protein_id": ["P1"]},
                index=pd.Index(["P1;S1;"], name="site_id"),
            ),
            site_clusters=pd.Series(
                [1],
                index=pd.Index(["P1;S1;"], name="site_id"),
                dtype="int64",
            ),
            site_metadata=pd.DataFrame(
                {"gene_symbol": ["P1"]},
                index=pd.Index(["P1;S1;"], name="site_id"),
            ),
            prediction_matrix=pd.DataFrame(
                {"K1": [0.9]},
                index=pd.Index(["P1;S1;"], name="site_id"),
            ),
            kinase_substrates={"K1": ("P1;S1;",)},
            substrate_support_cutoff=0.5,
            assignment_policy="cutoff_binary",
        )


def test_protein_context_builder_rejects_non_empty_missing_required_columns() -> None:
    with pytest.raises(WorkflowStageError, match="missing columns"):
        build_protein_site_context_table(
            site_membership=pd.DataFrame({"site_id": ["P1;S1;"]}),
        )

from __future__ import annotations

from dataclasses import fields
from pathlib import Path

import pandas as pd
import pytest

from phospy import (
    AnalysisReadyPhosphoDataset,
    KinaseWorkflow,
    SignalomeWorkflow,
)
from phospy.api import (
    KinasePredictionConfig,
    KinaseScoringConfig,
    KinaseWorkflowRequest,
    Organism,
    ReferenceBundle,
    ReferenceContextCompatibilityPolicy,
    SignalomeConfig,
    SignalomeWorkflowRequest,
)
from phospy.api.configs import SignalomeOutputConfig
from phospy.science.signalomes.science import (
    build_kinase_network,
    build_module_assignments,
)
from phospy.workflows.kinase.science import build_kinase_profiles
from tests.support.intensity_scale_states import (
    supported_linear_intensity_scale_state,
    supported_linear_processing_state,
)
from tests.support.rewrite_fixture_data import (
    load_adaptive_sampling_edge_rank_weighted_fusion_scores,
)
from tests.support.signalome_config import build_signalome_config
from tests.support.site_keys import (
    site_key_context_columns,
    site_key_index_from_display_ids,
)

ROOT = Path(__file__).resolve().parents[2]


def _dataset() -> AnalysisReadyPhosphoDataset:
    display_ids = ["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"]
    site_index = site_key_index_from_display_ids(display_ids)
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 4.0],
            "sample_b": [2.0, 4.0, 1.0],
        },
        index=site_index,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            **site_key_context_columns(site_index),
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": [
                ("A" * 15) + str(site).strip().upper()[0] + ("A" * 15)
                for site in ["Y182", "S9", "T308"]
            ],
            "protein_id": ["MAPK14", "GSK3B", "AKT1"],
        },
        index=site_index.copy(),
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


def _references() -> ReferenceBundle:
    site_ids = pd.Index(["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"], name="site_id")
    return ReferenceBundle(
        organism=Organism.RAT,
        kinase_substrate_map=pd.DataFrame(
            {
                "kinase": ["MAP2K6", "MAP2K6", "AKT1", "AKT1"],
                "substrate_site": [
                    "MAPK14;Y182;",
                    "GSK3B;S9;",
                    "GSK3B;S9;",
                    "AKT1;T308;",
                ],
            }
        ),
        site_sequences=pd.DataFrame(
            {"site_sequence": ["A" * 31, "B" * 31, "C" * 31]},
            index=site_ids,
        ),
    )


def test_profile_policy_historical_baseline_locks_strict_median_behavior_and_contract_surface() -> (
    None
):
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 3.0],
            "sample_b": [2.0, 4.0],
            "sample_c": [3.0, 5.0],
        },
        index=pd.Index(["SITE_1", "SITE_2"], name="site_id"),
    )
    kinase_substrate_map = pd.DataFrame(
        {
            "kinase": ["K1", "K1"],
            "substrate_site": ["SITE_1", "SITE_2"],
        }
    )

    result = build_kinase_profiles(
        phospho=phospho,
        kinase_substrate_map=kinase_substrate_map,
        min_substrates=2,
    )

    assert result.profile_matrix.at["K1", "sample_c"] == pytest.approx(4.0)
    phospho_with_missing = phospho.copy(deep=True)
    phospho_with_missing.loc["SITE_2", "sample_c"] = float("nan")
    strict = build_kinase_profiles(
        phospho=phospho_with_missing,
        kinase_substrate_map=kinase_substrate_map,
        min_substrates=2,
        profile_missing_value_strategy="strict",
    )
    median_skipna = build_kinase_profiles(
        phospho=phospho_with_missing,
        kinase_substrate_map=kinase_substrate_map,
        min_substrates=2,
        profile_missing_value_strategy="median_skipna",
    )
    assert pd.isna(strict.profile_matrix.at["K1", "sample_c"])
    assert median_skipna.profile_matrix.at["K1", "sample_c"] == pytest.approx(3.0)
    assert "profile_missing_value_strategy" in {
        field.name for field in fields(KinaseScoringConfig)
    }
    assert KinaseScoringConfig().profile_missing_value_strategy == "strict"
    with pytest.raises(TypeError, match="profile_policy"):
        KinaseScoringConfig(  # type: ignore[call-arg]
            min_substrates=2,
            profile_policy={"missing_value_strategy": "median_skipna"},
        )


def test_adaptive_sampling_historical_baseline_is_archival_and_svm_mode_is_not_rewrite_contract() -> (
    None
):
    prediction_fields = {field.name for field in fields(KinasePredictionConfig)}
    assert prediction_fields == {
        "top_k",
        "deterministic_max_selected_kinases",
        "adaptive_ensemble_runs",
        "mode",
        "adaptive_policy",
        "n_iterations",
        "random_state",
    }
    assert prediction_fields.isdisjoint(
        {
            "svm_mode",
            "allow_profile_only_fallback",
            "score_threshold",
            "inclusion",
            "min_motif_size",
        }
    )
    with pytest.raises(TypeError, match="svm_mode"):
        KinasePredictionConfig(  # type: ignore[call-arg]
            top_k=3,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
            svm_mode="r_parity",
        )
    with pytest.raises(TypeError, match="score_threshold"):
        KinasePredictionConfig(  # type: ignore[call-arg]
            top_k=3,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
            score_threshold=0.8,
        )
    with pytest.raises(TypeError, match="allow_profile_only_fallback"):
        KinasePredictionConfig(  # type: ignore[call-arg]
            top_k=3,
            deterministic_max_selected_kinases=2,
            adaptive_ensemble_runs=2,
            allow_profile_only_fallback=True,
        )
    assert KinasePredictionConfig().mode == "deterministic_ranking"

    promoted_trace_scores = load_adaptive_sampling_edge_rank_weighted_fusion_scores()
    assert not promoted_trace_scores.empty
    assert list(promoted_trace_scores.index.astype(str))[:2] == ["SITE_A", "SITE_B"]


def test_signalome_clustering_historical_baseline_locks_dominant_module_assignment_behavior() -> (
    None
):
    display_ids = ["P1;S1;", "P1;S2;", "P2;S3;", "P2;S4;"]
    site_index = site_key_index_from_display_ids(display_ids)
    prediction_matrix = pd.DataFrame(
        {
            "KINASE_A": [0.9, 0.1, 0.1, 0.2],
            "KINASE_B": [0.8, 0.95, 0.8, 0.85],
        },
        index=site_index.copy(),
        dtype=float,
    )
    site_to_protein = pd.Series(
        ["P1", "P1", "P2", "P2"],
        index=prediction_matrix.index.copy(),
        name="protein_id",
        dtype=str,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            "gene_symbol": ["P1", "P1", "P2", "P2"],
            "site": ["S1", "S2", "S3", "S4"],
        },
        index=site_index.copy(),
    )

    assignments = build_module_assignments(
        prediction_matrix=prediction_matrix,
        site_to_protein=site_to_protein,
        site_metadata=site_metadata,
    )
    p1s1, p1s2, p2s3, p2s4 = site_index.astype(str).tolist()

    assert assignments.loc[p1s1, "module_top_kinase"] == "KINASE_A"
    assert assignments.loc[p1s2, "module_top_kinase"] == "KINASE_A"
    assert int(assignments.loc[p1s1, "module_id"]) == 1
    assert int(assignments.loc[p1s2, "module_id"]) == 1
    assert assignments.loc[p2s3, "module_top_kinase"] == "KINASE_B"
    assert assignments.loc[p2s4, "module_top_kinase"] == "KINASE_B"
    assert int(assignments.loc[p2s3, "module_id"]) == 2
    with pytest.raises(TypeError, match="module_selection_policy"):
        SignalomeConfig(  # type: ignore[call-arg]
            module_selection_policy={"strategy": "single_module"},
        )


def test_weighted_top_assignment_historical_baseline_locks_fractional_metadata_and_non_fractional_module_selection() -> (
    None
):
    display_ids = ["P1;S1;", "P1;S2;"]
    site_index = site_key_index_from_display_ids(display_ids)
    prediction_matrix = pd.DataFrame(
        {
            "KINASE_A": [0.95, 0.2],
            "KINASE_B": [0.95, 0.96],
        },
        index=site_index.copy(),
        dtype=float,
    )
    site_to_protein = pd.Series(
        ["P1", "P1"],
        index=prediction_matrix.index.copy(),
        name="protein_id",
        dtype=str,
    )
    site_metadata = pd.DataFrame(
        {
            "site_key": site_index.astype(str).tolist(),
            "display_id": display_ids,
            "gene_symbol": ["P1", "P1"],
            "site": ["S1", "S2"],
        },
        index=site_index.copy(),
    )

    assignments = build_module_assignments(
        prediction_matrix=prediction_matrix,
        site_to_protein=site_to_protein,
        site_metadata=site_metadata,
    )
    p1s1, p1s2 = site_index.astype(str).tolist()
    tied = assignments.loc[p1s1]

    assert tied["top_kinase_candidates"] == ("KINASE_A", "KINASE_B")
    assert tied["top_kinase_weights"] == (("KINASE_A", 0.5), ("KINASE_B", 0.5))
    assert sum(weight for _, weight in tied["top_kinase_weights"]) == pytest.approx(1.0)
    assert assignments.loc[p1s1, "module_top_kinase"] == "KINASE_A"
    assert assignments.loc[p1s2, "module_top_kinase"] == "KINASE_A"
    assert int(assignments.loc[p1s1, "module_id"]) == 1


def test_network_policy_variant_historical_baseline_locks_signed_edges_and_narrow_config_surface() -> (
    None
):
    downstream_scores = pd.DataFrame(
        {
            "KINASE_A": [1.0, 2.0, 3.0, 4.0],
            "KINASE_B": [4.0, 3.0, 2.0, 1.0],
            "KINASE_C": [1.0, 1.0, 1.0, 1.0],
        },
        index=pd.Index(["S1", "S2", "S3", "S4"], name="site_id"),
    )
    positive_only_edges, _ = build_kinase_network(
        downstream_score_matrix=downstream_scores,
        kinase_order=["KINASE_A", "KINASE_B", "KINASE_C"],
        kinase_substrates={
            "KINASE_A": ("S1", "S2"),
            "KINASE_B": ("S3", "S4"),
            "KINASE_C": ("S1",),
        },
        threshold=0.9,
        network_policy="positive_only",
    )
    absolute_threshold_edges, _ = build_kinase_network(
        downstream_score_matrix=downstream_scores,
        kinase_order=["KINASE_A", "KINASE_B", "KINASE_C"],
        kinase_substrates={
            "KINASE_A": ("S1", "S2"),
            "KINASE_B": ("S3", "S4"),
            "KINASE_C": ("S1",),
        },
        threshold=0.9,
        network_policy="absolute_threshold",
    )
    signed_edges, _ = build_kinase_network(
        downstream_score_matrix=downstream_scores,
        kinase_order=["KINASE_A", "KINASE_B", "KINASE_C"],
        kinase_substrates={
            "KINASE_A": ("S1", "S2"),
            "KINASE_B": ("S3", "S4"),
            "KINASE_C": ("S1",),
        },
        threshold=0.9,
        network_policy="signed",
    )

    assert positive_only_edges.empty
    assert absolute_threshold_edges.shape[0] == 1
    assert absolute_threshold_edges.at[0, "source_kinase"] == "KINASE_A"
    assert absolute_threshold_edges.at[0, "target_kinase"] == "KINASE_B"
    assert absolute_threshold_edges.at[0, "correlation"] == pytest.approx(1.0)
    assert signed_edges.shape[0] == 1
    assert signed_edges.at[0, "source_kinase"] == "KINASE_A"
    assert signed_edges.at[0, "target_kinase"] == "KINASE_B"
    assert signed_edges.at[0, "correlation"] == pytest.approx(-1.0)

    assert "output" in {field.name for field in fields(SignalomeConfig)}
    assert SignalomeConfig().output.network_policy == "signed"
    assert SignalomeConfig(
        output=SignalomeOutputConfig(network_policy="positive_only")
    ).output.network_policy == ("positive_only")
    assert SignalomeConfig().scientific.assignment_policy == "cutoff_binary"
    assert SignalomeConfig().clustering.clustering_engine == "scipy_hierarchical"
    assert SignalomeConfig().performance.max_exact_tree_sites == 2000
    assert SignalomeConfig().validation.score_preconditioning_policy == (
        "error_on_drop"
    )


def test_expanded_signalome_historical_baseline_locks_supported_lane_to_materialized_output() -> (
    None
):
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
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
            ),
            activity_config=None,
        )
    )
    signalome_result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=build_signalome_config(
                substrate_support_cutoff=0.5,
                network_correlation_threshold=0.5,
                reference_context_compatibility_policy=(
                    ReferenceContextCompatibilityPolicy.ALLOW_UNKNOWN_WITH_CAVEAT
                ),
            ),
        )
    )

    expanded = signalome_result.expanded_signalome
    assert expanded is not None
    assert not expanded.empty
    assert {"kinase", "row_kind", "site_id", "site_order"}.issubset(
        set(expanded.columns)
    )


def test_activity_parity_lock_historical_baseline_uses_rewrite_owned_fixture_path() -> (
    None
):
    source = (ROOT / "tests" / "parity" / "test_activity_stage_parity.py").read_text(
        encoding="utf-8"
    )
    assert "from tests.support.rewrite_fixture_data import (" in source
    assert "load_activity_reference_predmat" in source
    assert "load_activity_reference_activity_matrix" in source
    assert "load_activity_reference_thresholded_substrate_mean_activity" in source
    assert "load_activity_reference_target_table" in source
    assert "legacy_archive" not in source
    assert "tests_legacy/fixtures" not in source
    assert "pytest.mark.parity" in source
    assert "pytest.mark.activity_parity" in source


def test_activity_parity_lock_historical_baseline_requires_explicit_ci_gate() -> None:
    ci_workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )
    assert "activity-parity-gate:" in ci_workflow
    assert "tests/parity/test_activity_stage_parity.py" in ci_workflow
    assert "activity_parity" in ci_workflow

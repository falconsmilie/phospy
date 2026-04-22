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
    SignalomeConfig,
    SignalomeWorkflowRequest,
)
from phospy.signalomes.science import (
    build_kinase_network,
    build_module_assignments,
)
from phospy.workflows.kinase.science import build_kinase_profiles
from tests.support.rewrite_fixture_data import (
    load_adaptive_sampling_edge_combined_scores,
)
from tests.support.transformation_states import supported_linear_state

ROOT = Path(__file__).resolve().parents[2]


def _dataset() -> AnalysisReadyPhosphoDataset:
    site_ids = pd.Index(["MAPK14;Y182;", "GSK3B;S9;", "AKT1;T308;"], name="site_id")
    phospho = pd.DataFrame(
        {
            "sample_a": [1.0, 2.0, 4.0],
            "sample_b": [2.0, 4.0, 1.0],
        },
        index=site_ids,
    )
    site_metadata = pd.DataFrame(
        {
            "gene_symbol": ["MAPK14", "GSK3B", "AKT1"],
            "site": ["Y182", "S9", "T308"],
            "site_sequence": ["A" * 31, "B" * 31, "C" * 31],
            "protein_id": ["MAPK14", "GSK3B", "AKT1"],
        },
        index=site_ids,
    )
    return AnalysisReadyPhosphoDataset(
        phospho=phospho,
        site_metadata=site_metadata,
        organism=Organism.RAT,
        transformation_state=supported_linear_state(has_total_matrix=False),
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


def test_profile_policy_donor_locks_strict_median_behavior_and_contract_surface() -> (
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


def test_adaptive_sampling_donor_is_archival_and_svm_mode_is_not_rewrite_contract() -> (
    None
):
    prediction_fields = {field.name for field in fields(KinasePredictionConfig)}
    assert prediction_fields == {
        "top_k",
        "ensemble_size",
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
            ensemble_size=2,
            svm_mode="r_parity",
        )
    with pytest.raises(TypeError, match="score_threshold"):
        KinasePredictionConfig(  # type: ignore[call-arg]
            top_k=3,
            ensemble_size=2,
            score_threshold=0.8,
        )
    with pytest.raises(TypeError, match="allow_profile_only_fallback"):
        KinasePredictionConfig(  # type: ignore[call-arg]
            top_k=3,
            ensemble_size=2,
            allow_profile_only_fallback=True,
        )
    assert KinasePredictionConfig().mode == "deterministic_ranking"

    promoted_trace_scores = load_adaptive_sampling_edge_combined_scores()
    assert not promoted_trace_scores.empty
    assert list(promoted_trace_scores.index.astype(str))[:2] == ["SITE_A", "SITE_B"]


def test_signalome_clustering_donor_locks_rewrite_dominant_module_assignment_behavior() -> (
    None
):
    prediction_matrix = pd.DataFrame(
        {
            "KINASE_A": [0.9, 0.1, 0.1, 0.2],
            "KINASE_B": [0.8, 0.95, 0.8, 0.85],
        },
        index=pd.Index(["P1;S1;", "P1;S2;", "P2;S3;", "P2;S4;"], name="site_id"),
        dtype=float,
    )
    site_to_protein = pd.Series(
        ["P1", "P1", "P2", "P2"],
        index=prediction_matrix.index.copy(),
        name="protein_id",
        dtype=str,
    )

    assignments = build_module_assignments(
        prediction_matrix=prediction_matrix,
        site_to_protein=site_to_protein,
    )

    assert assignments.loc["P1;S1;", "module_top_kinase"] == "KINASE_A"
    assert assignments.loc["P1;S2;", "module_top_kinase"] == "KINASE_A"
    assert int(assignments.loc["P1;S1;", "module_id"]) == 1
    assert int(assignments.loc["P1;S2;", "module_id"]) == 1
    assert assignments.loc["P2;S3;", "module_top_kinase"] == "KINASE_B"
    assert assignments.loc["P2;S4;", "module_top_kinase"] == "KINASE_B"
    assert int(assignments.loc["P2;S3;", "module_id"]) == 2
    with pytest.raises(TypeError, match="module_selection_policy"):
        SignalomeConfig(  # type: ignore[call-arg]
            substrate_support_cutoff=0.5,
            network_correlation_threshold=0.5,
            module_selection_policy={"strategy": "single_module"},
        )


def test_weighted_top_assignment_donor_locks_fractional_metadata_and_non_fractional_module_selection() -> (
    None
):
    prediction_matrix = pd.DataFrame(
        {
            "KINASE_A": [0.95, 0.2],
            "KINASE_B": [0.95, 0.96],
        },
        index=pd.Index(["P1;S1;", "P1;S2;"], name="site_id"),
        dtype=float,
    )
    site_to_protein = pd.Series(
        ["P1", "P1"],
        index=prediction_matrix.index.copy(),
        name="protein_id",
        dtype=str,
    )

    assignments = build_module_assignments(
        prediction_matrix=prediction_matrix,
        site_to_protein=site_to_protein,
    )
    tied = assignments.loc["P1;S1;"]

    assert tied["top_kinase_candidates"] == ("KINASE_A", "KINASE_B")
    assert tied["top_kinase_weights"] == (("KINASE_A", 0.5), ("KINASE_B", 0.5))
    assert sum(weight for _, weight in tied["top_kinase_weights"]) == pytest.approx(1.0)
    assert assignments.loc["P1;S1;", "module_top_kinase"] == "KINASE_A"
    assert assignments.loc["P1;S2;", "module_top_kinase"] == "KINASE_A"
    assert int(assignments.loc["P1;S1;", "module_id"]) == 1


def test_network_policy_variant_donor_locks_signed_edges_and_narrow_config_surface() -> (
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

    assert "network_policy" in {field.name for field in fields(SignalomeConfig)}
    assert SignalomeConfig().network_policy == "signed"
    assert SignalomeConfig(network_policy="positive_only").network_policy == (
        "positive_only"
    )


def test_expanded_signalome_donor_locks_supported_lane_to_materialized_output() -> None:
    kinase_result = KinaseWorkflow().run(
        KinaseWorkflowRequest(
            dataset=_dataset(),
            references=_references(),
            scoring_config=KinaseScoringConfig(min_substrates=2),
            prediction_config=KinasePredictionConfig(top_k=2, ensemble_size=2),
            activity_config=None,
        )
    )
    signalome_result = SignalomeWorkflow().run(
        SignalomeWorkflowRequest(
            kinase_result=kinase_result,
            config=SignalomeConfig(
                substrate_support_cutoff=0.5,
                network_correlation_threshold=0.5,
            ),
        )
    )

    expanded = signalome_result.expanded_signalome
    assert expanded is not None
    assert not expanded.empty
    assert {"kinase", "row_kind", "site_id", "site_order"}.issubset(
        set(expanded.columns)
    )


def test_activity_parity_lock_donor_uses_rewrite_owned_fixture_path() -> None:
    source = (ROOT / "tests" / "parity" / "test_activity_stage_parity.py").read_text(
        encoding="utf-8"
    )
    assert "from tests.support.rewrite_fixture_data import (" in source
    assert "load_activity_reference_predmat" in source
    assert "load_activity_reference_weighted_activity" in source
    assert "load_activity_reference_ksea_scores" in source
    assert "load_activity_reference_target_table" in source
    assert "legacy_archive" not in source
    assert "tests_legacy/fixtures" not in source
    assert "pytest.mark.parity" in source
    assert "pytest.mark.activity_parity" in source


def test_activity_parity_lock_donor_requires_explicit_ci_gate() -> None:
    ci_workflow = (ROOT / ".github" / "workflows" / "ci.yml").read_text(
        encoding="utf-8"
    )
    assert "activity-parity-gate:" in ci_workflow
    assert "tests/parity/test_activity_stage_parity.py" in ci_workflow
    assert "activity_parity" in ci_workflow

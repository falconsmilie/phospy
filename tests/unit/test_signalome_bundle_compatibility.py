from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import SignalomeConfig
from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
    SignalomeOutputConfig,
    SignalomeScientificConfig,
    SignalomeValidationConfig,
)
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._signalome.compatibility import (
    normalize_module_assignments_table,
    signalome_alignment_diagnostics_from_payload,
    signalome_alignment_diagnostics_to_payload,
    signalome_module_selection_diagnostics_from_payload,
    signalome_module_selection_diagnostics_to_payload,
    signalome_network_correlation_diagnostics_from_payload,
    signalome_network_correlation_diagnostics_to_payload,
    signalome_score_preconditioning_diagnostics_from_payload,
    signalome_score_preconditioning_diagnostics_to_payload,
)
from phospy.io.bundles.signalome import SignalomeWorkflowConfigSnapshot
from phospy.science.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    SignalomeAlignmentDiagnostics,
    SignalomeAlignmentInputDiagnostics,
    SignalomeModuleSelectionDiagnostics,
    SignalomeNetworkCorrelationDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
)


def _full_signalome_snapshot_payload(
    **overrides: object,
) -> dict[str, dict[str, object]]:
    signalome_config: dict[str, object] = {
        "scientific": {
            "substrate_support_cutoff": 0.5,
            "assignment_policy": SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
        },
        "clustering": {
            "module_count": None,
            "module_selection_primary_correlation_threshold": 0.6,
            "module_selection_fallback_correlation_threshold": 0.2,
            "module_selection_max_clusters": 15,
            "candidate_scoring_policy": "full",
            "clustering_engine": "exact_python",
        },
        "validation": {
            "score_preconditioning_policy": (
                SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
            ),
            "allow_mixed_total_protein_quantitative_meaning": False,
        },
        "output": {
            "network_correlation_threshold": 0.6,
            "network_policy": "signed",
        },
        "performance": {
            "max_exact_tree_sites": 2000,
            "max_full_candidate_scoring_sites": 2000,
        },
    }
    for key, value in overrides.items():
        if key in {"substrate_support_cutoff", "assignment_policy"}:
            scientific = signalome_config["scientific"]
            assert isinstance(scientific, dict)
            scientific[key] = value
            continue
        if key in {
            "module_count",
            "module_selection_primary_correlation_threshold",
            "module_selection_fallback_correlation_threshold",
            "module_selection_max_clusters",
            "candidate_scoring_policy",
            "clustering_engine",
        }:
            clustering = signalome_config["clustering"]
            assert isinstance(clustering, dict)
            clustering[key] = value
            continue
        if key in {"score_preconditioning_policy"}:
            validation = signalome_config["validation"]
            assert isinstance(validation, dict)
            validation[key] = value
            continue
        if key in {"allow_mixed_total_protein_quantitative_meaning"}:
            validation = signalome_config["validation"]
            assert isinstance(validation, dict)
            validation[key] = value
            continue
        if key in {"network_correlation_threshold", "network_policy"}:
            output = signalome_config["output"]
            assert isinstance(output, dict)
            output[key] = value
            continue
        if key in {"max_exact_tree_sites", "max_full_candidate_scoring_sites"}:
            performance = signalome_config["performance"]
            assert isinstance(performance, dict)
            performance[key] = value
            continue
        signalome_config[key] = value
    return {"signalome_config": signalome_config}


def test_signalome_snapshot_rejects_partial_payload_without_required_fields() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="config snapshot.signalome_config is missing required field\\(s\\):",
    ):
        SignalomeWorkflowConfigSnapshot.from_payload({"signalome_config": {}})


def test_signalome_snapshot_rejects_removed_signalome_cutoff_alias() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="contains unsupported field\\(s\\): signalome_cutoff",
    ):
        SignalomeWorkflowConfigSnapshot.from_payload(
            {"signalome_config": {"signalome_cutoff": 0.6}}
        )


def test_signalome_snapshot_supports_assignment_policy_payload() -> None:
    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(
        _full_signalome_snapshot_payload(
            assignment_policy=SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP
        )
    )

    assert (
        snapshot.signalome_config.scientific.assignment_policy
        == SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP
    )


def test_signalome_snapshot_supports_network_policy_payload() -> None:
    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(
        _full_signalome_snapshot_payload(network_policy="absolute_threshold")
    )

    assert snapshot.signalome_config.output.network_policy == "absolute_threshold"


def test_signalome_snapshot_round_trip_preserves_network_policy() -> None:
    snapshot = SignalomeWorkflowConfigSnapshot(
        signalome_config=SignalomeConfig(
            scientific=SignalomeScientificConfig(
                substrate_support_cutoff=0.5,
                assignment_policy=SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
            ),
            validation=SignalomeValidationConfig(
                score_preconditioning_policy=(
                    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
                ),
            ),
            output=SignalomeOutputConfig(
                network_correlation_threshold=0.7,
                network_policy="absolute_threshold",
            ),
        )
    )

    restored = SignalomeWorkflowConfigSnapshot.from_payload(snapshot.to_payload())
    assert restored == snapshot


def test_signalome_snapshot_payload_round_trip_preserves_all_fields() -> None:
    payload = {
        "signalome_config": {
            "scientific": {
                "substrate_support_cutoff": 0.42,
                "assignment_policy": SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
            },
            "clustering": {
                "module_count": 6,
                "module_selection_primary_correlation_threshold": 0.67,
                "module_selection_fallback_correlation_threshold": 0.23,
                "module_selection_max_clusters": 15,
                "candidate_scoring_policy": "sampled",
                "clustering_engine": "exact_python",
            },
            "validation": {
                "score_preconditioning_policy": (
                    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
                ),
                "allow_mixed_total_protein_quantitative_meaning": True,
            },
            "output": {
                "network_correlation_threshold": 0.73,
                "network_policy": "absolute_threshold",
            },
            "performance": {
                "max_exact_tree_sites": 2500,
                "max_full_candidate_scoring_sites": 1700,
            },
        }
    }

    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(payload)
    assert snapshot.to_payload() == payload


def test_signalome_snapshot_rejects_removed_max_exact_clustering_sites_alias() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="contains unsupported field\\(s\\): max_exact_clustering_sites",
    ):
        SignalomeWorkflowConfigSnapshot.from_payload(
            {
                "signalome_config": {
                    "max_exact_clustering_sites": 2500,
                }
            }
        )


def test_signalome_snapshot_rejects_unknown_clustering_engine_value() -> None:
    with pytest.raises(
        PhosPyInputError,
        match=(
            "config snapshot.signalome_config.clustering.clustering_engine "
            "must be one of:"
        ),
    ):
        SignalomeWorkflowConfigSnapshot.from_payload(
            _full_signalome_snapshot_payload(clustering_engine="approximate")
        )


def test_signalome_snapshot_accepts_engine_policy_fields() -> None:
    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(
        _full_signalome_snapshot_payload(
            candidate_scoring_policy="sampled",
            clustering_engine="scipy_hierarchical",
        )
    )

    assert snapshot.signalome_config.clustering.candidate_scoring_policy == "sampled"
    assert (
        snapshot.signalome_config.clustering.clustering_engine == "scipy_hierarchical"
    )


@pytest.mark.parametrize(
    "removed_name",
    [
        "cluster_tree_backend",
        "candidate_scoring_backend",
        "clustering_backend",
        "max_exact_cluster_tree_sites",
        "max_full_correlation_sites",
    ],
)
def test_signalome_snapshot_rejects_removed_backend_style_fields(
    removed_name: str,
) -> None:
    with pytest.raises(
        PhosPyInputError,
        match=f"contains unsupported field\\(s\\): {removed_name}",
    ):
        SignalomeWorkflowConfigSnapshot.from_payload(
            _full_signalome_snapshot_payload(
                **{removed_name: "full"},
            )
        )


def test_signalome_snapshot_rejects_float_module_count_payload() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="config snapshot.signalome_config.clustering.module_count must be an int",
    ):
        SignalomeWorkflowConfigSnapshot.from_payload(
            _full_signalome_snapshot_payload(module_count=6.0)
        )


def test_signalome_snapshot_rejects_removed_network_policy_alias() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="contains unsupported field\\(s\\): kinase_network_policy",
    ):
        SignalomeWorkflowConfigSnapshot.from_payload(
            {
                "signalome_config": {
                    "kinase_network_policy": "positive_only",
                }
            }
        )


def test_module_assignment_normalization_parses_serialized_fields() -> None:
    table = pd.DataFrame(
        {
            "module_id": [1, 2],
            "module_candidates": ["('AKT1', 'MAPK1')", ["GSK3B"]],
            "module_weights": ["{'AKT1': 0.4}", [("MAPK1", 0.9)]],
        }
    )

    normalized = normalize_module_assignments_table(table)

    assert list(normalized.loc[:, "module_candidates"]) == [
        ("AKT1", "MAPK1"),
        ("GSK3B",),
    ]
    assert list(normalized.loc[:, "module_weights"]) == [
        (("AKT1", 0.4),),
        (("MAPK1", 0.9),),
    ]


def test_module_selection_diagnostics_payload_round_trip() -> None:
    diagnostics = SignalomeModuleSelectionDiagnostics(
        strategy=SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
        selected_module_count=3,
        requested_module_count=None,
        threshold_used=0.5,
        max_clusters_evaluated=10,
        candidate_scores={},
        reason="selected primary threshold candidate",
        zero_variance_profile_count=1,
        near_constant_profile_count=2,
        excluded_from_correlation_count=3,
    )

    payload = signalome_module_selection_diagnostics_to_payload(diagnostics)
    restored = signalome_module_selection_diagnostics_from_payload(
        payload,
        scope="test",
    )

    assert restored == diagnostics


def test_score_preconditioning_diagnostics_payload_round_trip() -> None:
    diagnostics = SignalomeScorePreconditioningDiagnostics(
        policy=SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
        input_row_count=10,
        dropped_all_missing_row_count=2,
        retained_row_count=8,
    )

    payload = signalome_score_preconditioning_diagnostics_to_payload(diagnostics)
    restored = signalome_score_preconditioning_diagnostics_from_payload(
        payload,
        scope="test",
    )

    assert restored == diagnostics


def test_score_preconditioning_diagnostics_accepts_error_on_drop_policy() -> None:
    diagnostics = SignalomeScorePreconditioningDiagnostics(
        policy=SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
        input_row_count=10,
        dropped_all_missing_row_count=0,
        retained_row_count=10,
    )

    payload = signalome_score_preconditioning_diagnostics_to_payload(diagnostics)
    restored = signalome_score_preconditioning_diagnostics_from_payload(
        payload,
        scope="test",
    )

    assert restored == diagnostics


def test_score_preconditioning_diagnostics_requires_payload_mapping() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="test.score_preconditioning_diagnostics must be an object",
    ):
        signalome_score_preconditioning_diagnostics_from_payload(None, scope="test")


def test_network_correlation_diagnostics_payload_round_trip() -> None:
    diagnostics = SignalomeNetworkCorrelationDiagnostics(
        total_candidate_correlations=10,
        finite_correlations=3,
        undefined_correlations=7,
        constant_profile_correlations=2,
        insufficient_observation_correlations=1,
        missing_value_correlations=3,
        non_finite_value_correlations=1,
        edges_created=2,
        edges_skipped_non_finite_correlation=7,
    )

    payload = signalome_network_correlation_diagnostics_to_payload(diagnostics)
    restored = signalome_network_correlation_diagnostics_from_payload(
        payload,
        scope="test",
    )

    assert restored == diagnostics


def test_network_correlation_diagnostics_requires_payload_mapping() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="test.network_correlation_diagnostics must be an object",
    ):
        signalome_network_correlation_diagnostics_from_payload(None, scope="test")


def test_alignment_diagnostics_payload_round_trip() -> None:
    diagnostics = SignalomeAlignmentDiagnostics(
        dataset_sites=SignalomeAlignmentInputDiagnostics(
            provided_count=10,
            retained_count=8,
            dropped_count=2,
            dropped_reasons={
                "missing_from_prediction_scores": 1,
                "missing_from_downstream_scores": 1,
                "removed_by_score_preconditioning": 0,
                "removed_by_validation_policy": 0,
            },
        ),
        prediction_score_sites=SignalomeAlignmentInputDiagnostics(
            provided_count=11,
            retained_count=8,
            dropped_count=3,
            dropped_reasons={
                "missing_from_dataset": 2,
                "missing_from_downstream_scores": 1,
                "removed_by_score_preconditioning": 0,
                "removed_by_validation_policy": 0,
            },
        ),
        downstream_score_sites=SignalomeAlignmentInputDiagnostics(
            provided_count=9,
            retained_count=8,
            dropped_count=1,
            dropped_reasons={
                "missing_from_dataset": 0,
                "missing_from_prediction_scores": 0,
                "removed_by_score_preconditioning": 1,
                "removed_by_validation_policy": 0,
            },
        ),
        kinases=SignalomeAlignmentInputDiagnostics(
            provided_count=6,
            retained_count=4,
            dropped_count=2,
            dropped_reasons={
                "missing_from_prediction_scores": 1,
                "missing_from_downstream_scores": 1,
                "missing_kinase_support": 0,
            },
        ),
        protein_identifiers=SignalomeAlignmentInputDiagnostics(
            provided_count=8,
            retained_count=8,
            dropped_count=0,
            dropped_reasons={
                "removed_by_score_preconditioning": 0,
                "missing_protein_identifier": 0,
                "removed_by_validation_policy": 0,
            },
        ),
    )
    payload = signalome_alignment_diagnostics_to_payload(diagnostics)
    restored = signalome_alignment_diagnostics_from_payload(payload, scope="test")
    assert restored == diagnostics


def test_alignment_diagnostics_requires_payload_mapping() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="test.alignment_diagnostics must be an object",
    ):
        signalome_alignment_diagnostics_from_payload(None, scope="test")


def test_module_selection_diagnostics_rejects_partial_payload() -> None:
    payload = {
        "strategy": SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
        "selected_module_count": 3,
        "requested_module_count": None,
        "threshold_used": 0.5,
        "max_clusters_evaluated": 10,
        "candidate_scores": {},
        "reason": "selected primary threshold candidate",
        "zero_variance_profile_count": 1,
        "excluded_from_correlation_count": 3,
    }
    with pytest.raises(
        PhosPyInputError,
        match=(
            "test.module_selection_diagnostics is missing required field\\(s\\): "
            "near_constant_profile_count"
        ),
    ):
        signalome_module_selection_diagnostics_from_payload(payload, scope="test")


def test_network_correlation_diagnostics_rejects_partial_payload() -> None:
    payload = {
        "total_candidate_correlations": 10,
        "finite_correlations": 3,
        "undefined_correlations": 7,
        "constant_profile_correlations": 2,
        "insufficient_observation_correlations": 1,
        "missing_value_correlations": 3,
        "non_finite_value_correlations": 1,
        "edges_created": 2,
    }
    with pytest.raises(
        PhosPyInputError,
        match=(
            "test.network_correlation_diagnostics is missing required field\\(s\\): "
            "edges_skipped_non_finite_correlation"
        ),
    ):
        signalome_network_correlation_diagnostics_from_payload(payload, scope="test")

from __future__ import annotations

import pandas as pd
import pytest

from phospy.api import SignalomeConfig
from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
)
from phospy.errors.input import PhosPyInputError
from phospy.io.bundles._signalome.compatibility import (
    normalize_module_assignments_table,
    signalome_module_selection_diagnostics_from_payload,
    signalome_module_selection_diagnostics_to_payload,
    signalome_network_correlation_diagnostics_from_payload,
    signalome_network_correlation_diagnostics_to_payload,
    signalome_score_preconditioning_diagnostics_from_payload,
    signalome_score_preconditioning_diagnostics_to_payload,
)
from phospy.io.bundles.signalome import SignalomeWorkflowConfigSnapshot
from phospy.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    SignalomeModuleSelectionDiagnostics,
    SignalomeNetworkCorrelationDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
)


def _full_signalome_snapshot_payload(
    **overrides: object,
) -> dict[str, dict[str, object]]:
    signalome_config: dict[str, object] = {
        "substrate_support_cutoff": 0.5,
        "network_correlation_threshold": 0.6,
        "network_policy": "signed",
        "assignment_policy": SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
        "score_preconditioning_policy": (
            SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
        ),
        "cluster_tree_backend": "exact",
        "candidate_scoring_backend": "full",
        "max_exact_cluster_tree_sites": 2000,
        "max_full_correlation_sites": 2000,
        "module_count": None,
        "module_selection_primary_correlation_threshold": 0.6,
        "module_selection_fallback_correlation_threshold": 0.2,
        "module_selection_max_clusters": 15,
    }
    signalome_config.update(overrides)
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
        snapshot.signalome_config.assignment_policy
        == SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP
    )


def test_signalome_snapshot_supports_network_policy_payload() -> None:
    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(
        _full_signalome_snapshot_payload(network_policy="absolute_threshold")
    )

    assert snapshot.signalome_config.network_policy == "absolute_threshold"


def test_signalome_snapshot_round_trip_preserves_network_policy() -> None:
    snapshot = SignalomeWorkflowConfigSnapshot(
        signalome_config=SignalomeConfig(
            substrate_support_cutoff=0.5,
            network_correlation_threshold=0.7,
            network_policy="absolute_threshold",
            assignment_policy=SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
            score_preconditioning_policy=(
                SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
            ),
        )
    )

    restored = SignalomeWorkflowConfigSnapshot.from_payload(snapshot.to_payload())
    assert restored == snapshot


def test_signalome_snapshot_payload_round_trip_preserves_all_fields() -> None:
    payload = {
        "signalome_config": {
            "substrate_support_cutoff": 0.42,
            "network_correlation_threshold": 0.73,
            "network_policy": "absolute_threshold",
            "assignment_policy": SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
            "score_preconditioning_policy": (
                SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP
            ),
            "cluster_tree_backend": "exact",
            "candidate_scoring_backend": "sampled",
            "max_exact_cluster_tree_sites": 2500,
            "max_full_correlation_sites": 1700,
            "module_count": 6,
            "module_selection_primary_correlation_threshold": 0.67,
            "module_selection_fallback_correlation_threshold": 0.23,
            "module_selection_max_clusters": 15,
        }
    }

    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(payload)
    assert snapshot.to_payload() == payload


def test_signalome_snapshot_rejects_removed_backend_alias_payload_fields() -> None:
    with pytest.raises(
        PhosPyInputError,
        match=(
            "contains unsupported field\\(s\\): "
            "clustering_backend, max_exact_clustering_sites"
        ),
    ):
        SignalomeWorkflowConfigSnapshot.from_payload(
            {
                "signalome_config": {
                    "substrate_support_cutoff": 0.42,
                    "network_correlation_threshold": 0.73,
                    "clustering_backend": "approximate",
                    "max_exact_clustering_sites": 2500,
                }
            }
        )


def test_signalome_snapshot_rejects_removed_network_policy_alias() -> None:
    with pytest.raises(
        PhosPyInputError,
        match="contains unsupported field\\(s\\): kinase_network_policy",
    ):
        SignalomeWorkflowConfigSnapshot.from_payload(
            {
                "signalome_config": {
                    "substrate_support_cutoff": 0.5,
                    "network_correlation_threshold": 0.6,
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

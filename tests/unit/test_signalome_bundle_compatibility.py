from __future__ import annotations

import pandas as pd

from phospy.api import SignalomeConfig
from phospy.api.configs import (
    SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY,
    SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
)
from phospy.io.bundles._signalome.compatibility import (
    normalize_module_assignments_table,
    signalome_module_selection_diagnostics_from_payload_with_compatibility_support,
    signalome_module_selection_diagnostics_to_payload,
    signalome_score_preconditioning_diagnostics_from_payload_with_compatibility_support,
    signalome_score_preconditioning_diagnostics_to_payload,
)
from phospy.io.bundles.signalome import SignalomeWorkflowConfigSnapshot
from phospy.signalomes.models import (
    SIGNALOME_MODULE_SELECTION_STRATEGY_CORRELATION_THRESHOLDS,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT,
    SIGNALOME_SCORE_PRECONDITIONING_POLICY_ERROR_ON_DROP,
    SignalomeModuleSelectionDiagnostics,
    SignalomeScorePreconditioningDiagnostics,
    default_signalome_score_preconditioning_diagnostics,
)


def test_signalome_snapshot_supports_compatibility_cutoff_payload() -> None:
    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(
        {"signalome_config": {"signalome_cutoff": 0.6}}
    )

    assert snapshot.signalome_config.substrate_support_cutoff == 0.6
    assert snapshot.signalome_config.network_correlation_threshold == 0.6
    assert snapshot.signalome_config.network_policy == "signed"
    assert (
        snapshot.signalome_config.assignment_policy
        == SIGNALOME_ASSIGNMENT_POLICY_CUTOFF_BINARY
    )
    assert snapshot.signalome_config.score_preconditioning_policy == (
        SIGNALOME_SCORE_PRECONDITIONING_POLICY_ALLOW_AND_REPORT
    )


def test_signalome_snapshot_supports_assignment_policy_payload() -> None:
    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(
        {
            "signalome_config": {
                "substrate_support_cutoff": 0.5,
                "network_correlation_threshold": 0.6,
                "assignment_policy": SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP,
            }
        }
    )

    assert (
        snapshot.signalome_config.assignment_policy
        == SIGNALOME_ASSIGNMENT_POLICY_WEIGHTED_TOP
    )


def test_signalome_snapshot_supports_network_policy_payload_and_compatibility_alias() -> (
    None
):
    explicit = SignalomeWorkflowConfigSnapshot.from_payload(
        {
            "signalome_config": {
                "substrate_support_cutoff": 0.5,
                "network_correlation_threshold": 0.6,
                "network_policy": "absolute_threshold",
            }
        }
    )
    compatibility_alias = SignalomeWorkflowConfigSnapshot.from_payload(
        {
            "signalome_config": {
                "substrate_support_cutoff": 0.5,
                "network_correlation_threshold": 0.6,
                "kinase_network_policy": "positive_only",
            }
        }
    )

    assert explicit.signalome_config.network_policy == "absolute_threshold"
    assert compatibility_alias.signalome_config.network_policy == "positive_only"


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
            "module_count": 6,
            "module_selection_primary_correlation_threshold": 0.67,
            "module_selection_fallback_correlation_threshold": 0.23,
            "module_selection_max_clusters": 15,
        }
    }

    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(payload)
    assert snapshot.to_payload() == payload


def test_module_assignment_compat_normalization_parses_serialized_fields() -> None:
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
    restored = (
        signalome_module_selection_diagnostics_from_payload_with_compatibility_support(
            payload,
            scope="test",
        )
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
    restored = signalome_score_preconditioning_diagnostics_from_payload_with_compatibility_support(
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
    restored = signalome_score_preconditioning_diagnostics_from_payload_with_compatibility_support(
        payload,
        scope="test",
    )

    assert restored == diagnostics


def test_score_preconditioning_diagnostics_compatibility_payload_defaults() -> None:
    restored = signalome_score_preconditioning_diagnostics_from_payload_with_compatibility_support(
        None,
        scope="test",
    )
    assert restored == default_signalome_score_preconditioning_diagnostics()

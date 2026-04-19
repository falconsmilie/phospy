from __future__ import annotations

import pandas as pd

from phospy.io.bundles._signalome.compatibility import (
    normalize_module_assignments_table,
    signalome_module_selection_diagnostics_from_payload_with_legacy_support,
    signalome_module_selection_diagnostics_to_payload,
)
from phospy.io.bundles.signalome import SignalomeWorkflowConfigSnapshot
from phospy.signalomes.models import SignalomeModuleSelectionDiagnostics


def test_signalome_snapshot_supports_legacy_cutoff_payload() -> None:
    snapshot = SignalomeWorkflowConfigSnapshot.from_payload(
        {"signalome_config": {"signalome_cutoff": 0.6}}
    )

    assert snapshot.signalome_config.substrate_support_cutoff == 0.6
    assert snapshot.signalome_config.network_correlation_threshold == 0.6


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
        strategy="correlation_thresholds",
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
    restored = signalome_module_selection_diagnostics_from_payload_with_legacy_support(
        payload,
        scope="test",
    )

    assert restored == diagnostics

from __future__ import annotations

import pandas as pd

from phospy.io.bundles._signalome.compatibility import (
    normalize_module_assignments_table,
)
from phospy.io.bundles.signalome import SignalomeWorkflowConfigSnapshot


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
